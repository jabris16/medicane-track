#************************************************************************************************************
#************************************************************************************************************
#************************************************************************************************************

# N.B. this code is inteded to be run over a singular NetCDF file, for code to use algorithm over multiple files please contact the author

#************************************************************************************************************

import os, sys
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib

#************************************************************************************************************

#name of ncdf file (fill with appropriate name)
filename = []

# detection criteria constants (fill in with appropriate values) 
pressure_treshold = []    # pressure gradient of feature point vs. all surrouding cells
wind_threshold = []    # average wind-speed value used as minimum threshold for detection
precip_threhsold = []   # average precipitation value used minimum threshold for detection

#************************************************************************************************************
# read in file and establish variables

nc = Dataset(filename)
mslp = nc.variables['x'][:]    # MSLP variable where x = MSLP header name in filename ncdf 
bg_mslp = np.mean(mslp)     # background MSLP variable
anomaly = np.subtract(mslp,bg_mslp)    # MSLP anomaly variable
precip = nc.variables['a'][:]    # precipitation variable where a = precipitation header name in filename ncdf 
wind = nc.variables['b'][:]    # wind-speed variable where b = wind-speed header name in filename ncdf 
time = nc.variables['c'][:]    # time variable where c = time header name in filename ncdf 
lat = nc.variables['y'][:]    # latitude variable where y = latitude header name in filename ncdf 
lon = nc.variables['z'][:]    # longitude varibale where z = longitude header name in filename ncdf 

#************************************************************************************************************
# creation of shapes,arrays and lists

anomaly_shape = np.shape(anomaly)
ident = np.zeros(anomaly_shape,dtype = 'i')    # empty array to later store identified feature points
pres_adj_total = np.zeros(anomaly_shape)    # empty array to later store feature point adjacent MSLP  values
pres_adj_av = np.zeros(anomaly_shape)    # empty array to later store feature point adjacent average MSLP values
precip_adj_total = np.zeros(anomaly_shape)    # empty array to later store feature point adjacent precipitation  values
precip_adj_av = np.zeros(anomaly_shape)    # empty array to later store feature point adjacent average precipitation values
wind_adj_total = np.zeros(anomaly_shape)    # empty array to later store feature point adjacent wind  values
wind_adj_av = np.zeros(anomaly_shape)    # empty array to later store feature point adjacent average wind values

# creation of empty lists where values will later be stored mid-loop to be added to final feature point dataframe

lat_list = []
lon_list = []
time_list = []
mslp_list = []
precip_list = []
wind_list = []

#************************************************************************************************************
# loop through all data points across space and time to establish medicane feature points and store in dataframe

for i in xrange(anomaly_shape[0]):
    for j in xrange(1,anomaly_shape[1]-1):    # first and last points excluded to account for boundaries of domain
        for k in xrange(1,anomaly_shape[2]-1):    # first and last points exluded to account for boundaries of domain
            # to assess negative gradient feature points only and that feature pointvs adjacent cells gradient exceeds defined pressure gradient threhsold
            if anomaly[i,j,k] < 0 and \  
            anomaly[i,j+1,k] - anomaly[i,j,k] > pressure_threshold and \
            anomaly[i,j,k+1] - anomaly[i,j,k] > pressure_threshold and \
            anomaly[i,j-1,k] - anomaly[i,j,k] > pressure_threshold and \
            anomaly[i,j,k-1] - anomaly[i,j,k] > pressure_threshold:
                ident[i,j,k] = 1
            else:
                pass
            
            if ident[i,j,k] == 1:     # storing feature point adjacent precipitation and wind-speed values
                precip_adj_total[i,j,k] = precip[i,j+1,k]
            + precip[i,j,k+1]
            + precip[i,j-1,k]
            + precip[i,j,k-1]
            precip_adj_av[i,j,k] = precip_adj_total[i,j,k]/4
                wind_adj_total[i,j,k] = wind[i,j+1,k]
            + wind[i,j,k+1]
            + wind[i,j-1,k]
            + wind[i,j,k-1]
            wind_adj_av[i,j,k] = wind_adj_total[i,j,k]/4
            
            if precip_adj_av[i,j,k] > 40 and \
            wind_adj_av[i,j,k] > 5:    # identified feature point characteristics which exceed defined thresholds stored in lists before being converted to numpy arrays
                lat_list.append(lat[j])
                lat_store = np.array(lat_list)
                lon_list.append(lon[k])
                lon_store = np.array(lon_list)
                time_list.append(time[i])
                time_store = np.array(time_list)
                mslp_list.append(anomaly[i,j,k])
                mslp_store = np.array(mslp_list)
                precip_list.append(precip[i,j,k])
                precip_store = np.array(precip_list)
                wind_list.append(wind[i,j,k])
                wind_store = np.array(wind_list)
                
                # create dataframe with all stored feature point information:

                feature_dataframe = pd.DataFrame({'Latitude':lat_store[:],'Longitude':lon_store[:],'Time step (days since 1/12/1850)':time_store[:],'MSLP anomaly (Pa)':mslp_store[:],'Precipitation (mm/day)':precip_store[:],'Wind-speed (m/sec)':wind_store[:]})
                
