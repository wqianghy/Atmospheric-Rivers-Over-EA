# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 10:09:38 2021
@author: Qiang Wang
"""

import netCDF4 as nc
import numpy as np
import pymannkendall as mk

#%% Create average / std /cv
ivt_average=np.zeros((71,241,321))
ivt_std    =np.zeros((71,241,321))
ivt_cv     =np.zeros((71,241,321))

iwv_average=np.zeros((71,241,321))
iwv_std    =np.zeros((71,241,321))
iwv_cv     =np.zeros((71,241,321))

for year_i in range(1950,2021,1):
    
    time_i=year_i-1950
    ivt_path='D:/AR_EXP/EXP0/IVT/IVT'+str(year_i)+'.nc'
    ivt_yeari=nc.Dataset(ivt_path)
    ivt_i_data=ivt_yeari.variables['IVT_YEAR']
    ivt_average[time_i,:,:]=np.mean(ivt_i_data,axis=0)
    ivt_std[time_i,:,:]=np.std(ivt_i_data,axis=0)
    ivt_cv[time_i,:,:]=ivt_std[time_i,:,:]/ivt_average[time_i,:,:]
    
    iwv_path='D:/AR_EXP/EXP0/IWV/IWV'+str(year_i)+'.nc'
    iwv_yeari=nc.Dataset(iwv_path);
    iwv_i_data=iwv_yeari.variables['IWV_YEAR'];
    iwv_average[time_i,:,:]=np.mean(iwv_i_data,axis=0)
    iwv_std[time_i,:,:]=np.std(iwv_i_data,axis=0);
    iwv_cv[time_i,:,:] =iwv_std[time_i,:,:]/iwv_average[time_i,:,:]
    
    print(year_i)
    
#%% output NC

ar_lat=ivt_yeari.variables['Ilat'];
ar_lon=ivt_yeari.variables['Ilon'];
ar_time=range(1950,2021,1);
outputIVT_IWV_path = 'D:/AR_EXP/EXP0/statis/IVT_IWV_ave_std_mk.nc';
outputIVT_IWV = nc.Dataset(outputIVT_IWV_path, mode="w", format='NETCDF4')

# Add Dimensions
time = outputIVT_IWV.createDimension('time', None)
lat = outputIVT_IWV.createDimension('lat', len(ar_lat))
lon = outputIVT_IWV.createDimension('lon', len(ar_lon))

# Add NetCDF Variables
times = outputIVT_IWV.createVariable('time', 'i4', ('time',))
lats = outputIVT_IWV.createVariable('lat', 'f8', ('lat',))
lons = outputIVT_IWV.createVariable('lon', 'f8', ('lon',))
IVT_AVE_NC = outputIVT_IWV.createVariable('IVT_AVE_NC', 'f8', ('time', 'lat', 'lon',))
IVT_STD_NC = outputIVT_IWV.createVariable('IVT_STD_NC', 'f8', ('time', 'lat', 'lon',))
IVT_CV_NC = outputIVT_IWV.createVariable('IVT_CV_NC', 'f8', ('time', 'lat', 'lon',))

IWV_AVE_NC = outputIVT_IWV.createVariable('IWV_AVE_NC', 'f8', ('time', 'lat', 'lon',))
IWV_STD_NC = outputIVT_IWV.createVariable('IWV_STD_NC', 'f8', ('time', 'lat', 'lon',))
IWV_CV_NC = outputIVT_IWV.createVariable('IWV_CV_NC', 'f8', ('time', 'lat', 'lon',))

# Assign Latitude and Longitude Values
times[:]=np.array(ar_time);
lats[:]=np.array(ar_lat);
lons[:]=np.array(ar_lon);

# Assign NetCDF Data Values
IVT_AVE_NC[:,:,:] = ivt_average
IVT_STD_NC[:,:,:] = ivt_std
IVT_CV_NC[:,:,:]  = ivt_cv
IWV_AVE_NC[:,:,:] = iwv_average
IWV_STD_NC[:,:,:] = iwv_std
IWV_CV_NC[:,:,:]  = iwv_cv

outputIVT_IWV.close()

#%% read average / std
ivt_iwv_ave_std_mk=nc.Dataset('D:/AR_EXP/EXP0/statis/IVT_IWV_ave_std_mk.nc')

ivt_ave_nc=np.array(ivt_iwv_ave_std_mk.variables['IVT_AVE_NC']);
ivt_std_nc=np.array(ivt_iwv_ave_std_mk.variables['IVT_STD_NC']);
ivt_cv_nc=np.array(ivt_iwv_ave_std_mk.variables['IVT_CV_NC']);

iwv_ave_nc=np.array(ivt_iwv_ave_std_mk.variables['IWV_AVE_NC']);
iwv_std_nc=np.array(ivt_iwv_ave_std_mk.variables['IWV_STD_NC']);
iwv_cv_nc=np.array(ivt_iwv_ave_std_mk.variables['IWV_CV_NC']);

ar_lat=ivt_iwv_ave_std_mk.variables['lat'];
ar_lon=ivt_iwv_ave_std_mk.variables['lon'];
ar_time=ivt_iwv_ave_std_mk.variables['time'];

#%% MK

mk_results=np.zeros((12,len(ar_lat),len(ar_lon))) 
# 0-1 ivt_ave z/ slope;2-3 ivt_std z/ slope;4-5 iwv_ave z/ slope;6-7 iwv_std z/ slope
for lat_i in range(len(ar_lat)):
    for lon_i in range(len(ar_lon)):
        
        ivt_ave_mk=mk.pre_whitening_modification_test(ivt_ave_nc[:,lat_i,lon_i]);
        ivt_std_mk=mk.pre_whitening_modification_test(ivt_std_nc[:,lat_i,lon_i]);
        ivt_cv_mk=mk.pre_whitening_modification_test(ivt_cv_nc[:,lat_i,lon_i]);
        
        iwv_ave_mk=mk.pre_whitening_modification_test(iwv_ave_nc[:,lat_i,lon_i]);
        iwv_std_mk=mk.pre_whitening_modification_test(iwv_std_nc[:,lat_i,lon_i]);
        iwv_cv_mk=mk.pre_whitening_modification_test(iwv_cv_nc[:,lat_i,lon_i]);
        
        mk_results[0,lat_i,lon_i]=ivt_ave_mk.z;
        mk_results[1,lat_i,lon_i]=ivt_ave_mk.slope;
        mk_results[2,lat_i,lon_i]=ivt_std_mk.z;
        mk_results[3,lat_i,lon_i]=ivt_std_mk.slope;
        mk_results[4,lat_i,lon_i]=ivt_cv_mk.z;
        mk_results[5,lat_i,lon_i]=ivt_cv_mk.slope;

        mk_results[6,lat_i,lon_i]=iwv_ave_mk.z;
        mk_results[7,lat_i,lon_i]=iwv_ave_mk.slope;
        mk_results[8,lat_i,lon_i]=iwv_std_mk.z;
        mk_results[9,lat_i,lon_i]=iwv_std_mk.slope;
        mk_results[10,lat_i,lon_i]=iwv_cv_mk.z;
        mk_results[11,lat_i,lon_i]=iwv_cv_mk.slope;

        print(lat_i,lon_i)


#%% Create  EXP 1_1 / EXP 2_1 / EXP 3_1
# Exp1_1 : remove trends of IVT/  IWV mean
# Exp2_1 : remove trends of IVT/  IWV cv
# Exp3_1 : fix k mean std

ivt_ave_trends=mk_results[1,:,:]
ivt_std_trends=mk_results[3,:,:]
ivt_cv_trends =mk_results[5,:,:]

iwv_ave_trends=mk_results[7,:,:]
iwv_std_trends=mk_results[9,:,:]
iwv_cv_trends =mk_results[11,:,:]

for year_i in range(1950,2021,1):
    
    time_i=year_i-1950;
    
    ivt_exp0_path='D:/AR_EXP/EXP0/IVT/IVT'+str(year_i)+'.nc';
    iwv_exp0_path='D:/AR_EXP/EXP0/IWV/IWV'+str(year_i)+'.nc';
    
    ivt_yeari=nc.Dataset(ivt_exp0_path);
    iwv_yeari=nc.Dataset(iwv_exp0_path);
    
    ivt_time=ivt_yeari.variables['Itime']
    ivt_lat =ivt_yeari.variables['Ilat']
    ivt_lon =ivt_yeari.variables['Ilon']
    
    ivt_i_data=np.array(ivt_yeari.variables['IVT_YEAR']);
    iwv_i_data=np.array(iwv_yeari.variables['IWV_YEAR']);
    
    # mean
    ivt_ave_exp0_yeari=np.mean(ivt_i_data,axis=0)
    iwv_ave_exp0_yeari=np.mean(iwv_i_data,axis=0)
    
    # std
    ivt_std_exp0_yeari=np.std(ivt_i_data,axis=0)
    iwv_std_exp0_yeari=np.std(iwv_i_data,axis=0)

    # cv
    ivt_cv_exp0_yeari=ivt_std_exp0_yeari/ivt_ave_exp0_yeari
    iwv_cv_exp0_yeari=iwv_std_exp0_yeari/iwv_ave_exp0_yeari

    # k - new : K=(ivt-ivt_ave)/(ivt_ave*cv)
    exp_k_ivt=(ivt_i_data-ivt_ave_exp0_yeari)/(ivt_ave_exp0_yeari*ivt_cv_exp0_yeari)
    exp_k_iwv=(iwv_i_data-iwv_ave_exp0_yeari)/(iwv_ave_exp0_yeari*iwv_cv_exp0_yeari);
    
    # 
    trends_ave_ivt=time_i*ivt_ave_trends
    trends_ave_iwv=time_i*iwv_ave_trends
    
    #  
    trends_cv_ivt=time_i*ivt_cv_trends
    trends_cv_iwv=time_i*iwv_cv_trends
    
    # exp1
    exp1_ivt_i_data = (ivt_ave_exp0_yeari-trends_ave_ivt) * (1 + exp_k_ivt*ivt_cv_exp0_yeari)
    exp1_iwv_i_data = (iwv_ave_exp0_yeari-trends_ave_iwv) * (1 + exp_k_iwv*iwv_cv_exp0_yeari)
    
    # exp2
    exp2_ivt_i_data = ivt_ave_exp0_yeari*(1+exp_k_ivt*(ivt_cv_exp0_yeari-trends_cv_ivt))
    exp2_iwv_i_data = iwv_ave_exp0_yeari*(1+exp_k_iwv*(iwv_cv_exp0_yeari-trends_cv_iwv))

    # exp3 fix k,remove both mean and std trends
    exp3_ivt_i_data = (ivt_ave_exp0_yeari-trends_ave_ivt) * (1 + exp_k_ivt*(ivt_cv_exp0_yeari-trends_cv_ivt))
    exp3_iwv_i_data = (iwv_ave_exp0_yeari-trends_ave_iwv) * (1 + exp_k_iwv*(iwv_cv_exp0_yeari-trends_cv_iwv))
    
    # output
    output_exp1_ivt_path = 'D:/AR_EXP/EXP1_1/IVT/IVT'+str(year_i)+'.nc';
    output_exp1_iwv_path = 'D:/AR_EXP/EXP1_1/IWV/IWV'+str(year_i)+'.nc';
    output_exp2_ivt_path = 'D:/AR_EXP/EXP2_1/IVT/IVT'+str(year_i)+'.nc';
    output_exp2_iwv_path = 'D:/AR_EXP/EXP2_1/IWV/IWV'+str(year_i)+'.nc';
    output_exp3_ivt_path = 'D:/AR_EXP/EXP3_1/IVT/IVT'+str(year_i)+'.nc';
    output_exp3_iwv_path = 'D:/AR_EXP/EXP3_1/IWV/IWV'+str(year_i)+'.nc';
    
    output_exp1_ivt = nc.Dataset(output_exp1_ivt_path, mode="w", format='NETCDF4')
    output_exp1_iwv = nc.Dataset(output_exp1_iwv_path, mode="w", format='NETCDF4')
    output_exp2_ivt = nc.Dataset(output_exp2_ivt_path, mode="w", format='NETCDF4')
    output_exp2_iwv = nc.Dataset(output_exp2_iwv_path, mode="w", format='NETCDF4')
    output_exp3_ivt = nc.Dataset(output_exp3_ivt_path, mode="w", format='NETCDF4')
    output_exp3_iwv = nc.Dataset(output_exp3_iwv_path, mode="w", format='NETCDF4')
    
    # output_exp1_ivt
    time = output_exp1_ivt.createDimension('time', None)
    lat  = output_exp1_ivt.createDimension('lat', len(ar_lat))
    lon  = output_exp1_ivt.createDimension('lon', len(ar_lon))
    times = output_exp1_ivt.createVariable('time', 'i4', ('time',))
    lats = output_exp1_ivt.createVariable('lat', 'f8', ('lat',))
    lons = output_exp1_ivt.createVariable('lon', 'f8', ('lon',))
    ivt  = output_exp1_ivt.createVariable('ivt', 'f8', ('time', 'lat', 'lon',))
    times[:]=np.array(ivt_time);
    lats[:]=np.array(ivt_lat);
    lons[:]=np.array(ivt_lon);
    ivt[:,:,:] = exp1_ivt_i_data
    output_exp1_ivt.close()
    
    # output_exp1_iwv
    time = output_exp1_iwv.createDimension('time', None)
    lat  = output_exp1_iwv.createDimension('lat', len(ar_lat))
    lon  = output_exp1_iwv.createDimension('lon', len(ar_lon))
    times = output_exp1_iwv.createVariable('time', 'i4', ('time',))
    lats = output_exp1_iwv.createVariable('lat', 'f8', ('lat',))
    lons = output_exp1_iwv.createVariable('lon', 'f8', ('lon',))
    iwv  = output_exp1_iwv.createVariable('iwv', 'f8', ('time', 'lat', 'lon',))
    times[:]=np.array(ivt_time);
    lats[:]=np.array(ivt_lat);
    lons[:]=np.array(ivt_lon);
    iwv[:,:,:] = exp1_iwv_i_data
    output_exp1_iwv.close()
    
    # output_exp2_ivt
    time = output_exp2_ivt.createDimension('time', None)
    lat  = output_exp2_ivt.createDimension('lat', len(ar_lat))
    lon  = output_exp2_ivt.createDimension('lon', len(ar_lon))
    times = output_exp2_ivt.createVariable('time', 'i4', ('time',))
    lats = output_exp2_ivt.createVariable('lat', 'f8', ('lat',))
    lons = output_exp2_ivt.createVariable('lon', 'f8', ('lon',))
    ivt  = output_exp2_ivt.createVariable('ivt', 'f8', ('time', 'lat', 'lon',))
    times[:]=np.array(ivt_time);
    lats[:]=np.array(ivt_lat);
    lons[:]=np.array(ivt_lon);
    ivt[:,:,:] = exp2_ivt_i_data
    output_exp2_ivt.close()

    # output_exp2_iwv
    time = output_exp2_iwv.createDimension('time', None)
    lat  = output_exp2_iwv.createDimension('lat', len(ar_lat))
    lon  = output_exp2_iwv.createDimension('lon', len(ar_lon))
    times = output_exp2_iwv.createVariable('time', 'i4', ('time',))
    lats = output_exp2_iwv.createVariable('lat', 'f8', ('lat',))
    lons = output_exp2_iwv.createVariable('lon', 'f8', ('lon',))
    iwv  = output_exp2_iwv.createVariable('iwv', 'f8', ('time', 'lat', 'lon',))
    times[:]=np.array(ivt_time);
    lats[:]=np.array(ivt_lat);
    lons[:]=np.array(ivt_lon);
    iwv[:,:,:] = exp2_iwv_i_data
    output_exp2_iwv.close()

    # output_exp3_ivt
    time = output_exp3_ivt.createDimension('time', None)
    lat  = output_exp3_ivt.createDimension('lat', len(ar_lat))
    lon  = output_exp3_ivt.createDimension('lon', len(ar_lon))
    times = output_exp3_ivt.createVariable('time', 'i4', ('time',))
    lats = output_exp3_ivt.createVariable('lat', 'f8', ('lat',))
    lons = output_exp3_ivt.createVariable('lon', 'f8', ('lon',))
    ivt  = output_exp3_ivt.createVariable('ivt', 'f8', ('time', 'lat', 'lon',))
    times[:]=np.array(ivt_time);
    lats[:]=np.array(ivt_lat);
    lons[:]=np.array(ivt_lon);
    ivt[:,:,:] = exp3_ivt_i_data
    output_exp3_ivt.close()
    
    # output_exp3_iwv
    time = output_exp3_iwv.createDimension('time', None)
    lat  = output_exp3_iwv.createDimension('lat', len(ar_lat))
    lon  = output_exp3_iwv.createDimension('lon', len(ar_lon))
    times = output_exp3_iwv.createVariable('time', 'i4', ('time',))
    lats = output_exp3_iwv.createVariable('lat', 'f8', ('lat',))
    lons = output_exp3_iwv.createVariable('lon', 'f8', ('lon',))
    iwv  = output_exp3_iwv.createVariable('iwv', 'f8', ('time', 'lat', 'lon',))
    times[:]=np.array(ivt_time);
    lats[:]=np.array(ivt_lat);
    lons[:]=np.array(ivt_lon);
    iwv[:,:,:] = exp3_iwv_i_data
    output_exp3_iwv.close()
 
    print(year_i)
