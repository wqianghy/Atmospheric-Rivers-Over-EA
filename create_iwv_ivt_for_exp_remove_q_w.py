# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 11:13:08 2021

@author: wqiang
"""

import netCDF4 as nc
import numpy as np
import pymannkendall as mk
from scipy.integrate import simps

#%% average

q_ave=np.zeros((71,15,241,321))
u_ave=np.zeros((71,15,241,321))
v_ave=np.zeros((71,15,241,321))

for year_i in range(1950,2021,1):
    
    q_yeari=nc.Dataset('I:/vairable_for_ar/humid/h'+str(year_i)+'.nc')
    u_yeari=nc.Dataset('I:/vairable_for_ar/uwind/uwind'+str(year_i)+'.nc')
    v_yeari=nc.Dataset('I:/vairable_for_ar/vwind/vwind'+str(year_i)+'.nc')
    
    for lev in range(15):
        
        q_ave[year_i-1950,lev,:,:]=np.mean(q_yeari.variables['q'][:,lev,:,240:],axis=(0))
        u_ave[year_i-1950,lev,:,:]=np.mean(u_yeari.variables['u'][:,lev,:,240:],axis=(0))
        v_ave[year_i-1950,lev,:,:]=np.mean(v_yeari.variables['v'][:,lev,:,240:],axis=(0))
        
        print(year_i,lev)
        
#%% mk of trends

mk_results=np.zeros((6,15,241,321)) 

for lat_i in range(241):
    for lon_i in range(321):
        for lev in range(15):
            
            q_ave_mk=mk.pre_whitening_modification_test(q_ave[:,lev,lat_i,lon_i])
            u_ave_mk=mk.pre_whitening_modification_test(u_ave[:,lev,lat_i,lon_i])
            v_ave_mk=mk.pre_whitening_modification_test(v_ave[:,lev,lat_i,lon_i])
            mk_results[0,lev,lat_i,lon_i]=q_ave_mk.z
            mk_results[1,lev,lat_i,lon_i]=q_ave_mk.slope
            mk_results[2,lev,lat_i,lon_i]=u_ave_mk.z
            mk_results[3,lev,lat_i,lon_i]=u_ave_mk.slope
            mk_results[4,lev,lat_i,lon_i]=v_ave_mk.z
            mk_results[5,lev,lat_i,lon_i]=v_ave_mk.slope

        print(lat_i,lon_i)
        
#%% function: remove the trends of q

def nc_minus_slope_q(year_i,q_trends,level,lat,lon):
    
    time_i=year_i-1950
    q_path='I:/vairable_for_ar/humid/h'+str(year_i)+'.nc'
    q_yeari=nc.Dataset(q_path)
    time_q = q_yeari.variables['time']
    output_q_path = 'G:/vairable_for_ar/humid_no_trends/h'+str(year_i)+'.nc'
    output_q = nc.Dataset(output_q_path, mode="w", format='NETCDF4')
    output_q.createDimension('q_time', None)
    output_q.createDimension('q_lev', len(level))
    output_q.createDimension('q_lat', len(lat))
    output_q.createDimension('q_lon', len(lon))
    q_times= output_q.createVariable('q_time', 'i4', ('q_time',))
    q_levs = output_q.createVariable('q_lev', 'i4', ('q_lev',))
    q_lats = output_q.createVariable('q_lat', 'f4', ('q_lat',))
    q_lons = output_q.createVariable('q_lon', 'f4', ('q_lon',))
    q      = output_q.createVariable('q', 'f4', ('q_time','q_lev', 'q_lat', 'q_lon',))
    q_times[:]=np.array(time_q)
    q_levs[:] =np.array(level)
    q_lats[:] =np.array(lat)
    q_lons[:] =np.array(lon)
    humid_time_lev=q_yeari.variables['q'][:,:,:,240:]
    q_time_lev_trends=time_i*q_trends
    q[:,:,:,:] = humid_time_lev-q_time_lev_trends
    
    q_yeari.close() 
    output_q.close() 

#%% function: remove the trends of u

def nc_minus_slope_u(year_i,u_trends,level,lat,lon):
    
    time_i=year_i-1950
    u_path='I:/vairable_for_ar/uwind/uwind'+str(year_i)+'.nc'
    u_yeari=nc.Dataset(u_path)
    time_u = u_yeari.variables['time']
    output_u_path = 'G:/vairable_for_ar/uwind_no_trends/u'+str(year_i)+'.nc'
    output_u = nc.Dataset(output_u_path, mode="w", format='NETCDF4')
    output_u.createDimension('u_time', None)
    output_u.createDimension('u_lev', len(level))
    output_u.createDimension('u_lat', len(lat))
    output_u.createDimension('u_lon', len(lon))
    u_times= output_u.createVariable('u_time', 'i4', ('u_time',))
    u_levs = output_u.createVariable('u_lev',  'i4', ('u_lev',))
    u_lats = output_u.createVariable('u_lat',  'f4', ('u_lat',))
    u_lons = output_u.createVariable('u_lon',  'f4', ('u_lon',))
    u      = output_u.createVariable('u', 'f4', ('u_time','u_lev', 'u_lat', 'u_lon',))
    u_times[:]=np.array(time_u)
    u_levs[:] =np.array(level)
    u_lats[:] =np.array(lat)
    u_lons[:] =np.array(lon)
    uwind_time_lev=u_yeari.variables['u'][:,:,:,240:]
    u_time_lev_trends=time_i*u_trends
    u[:,:,:,:] = uwind_time_lev-u_time_lev_trends
    
    u_yeari.close()
    output_u.close()
    
#%% function: remove the trends of v

def nc_minus_slope_v(year_i,v_trends,level,lat,lon):
    
    time_i=year_i-1950
    v_path='I:/vairable_for_ar/vwind/vwind'+str(year_i)+'.nc'
    v_yeari=nc.Dataset(v_path)
    time_v = v_yeari.variables['time']
    output_v_path = 'G:/vairable_for_ar/vwind_no_trends/v'+str(year_i)+'.nc'
    output_v = nc.Dataset(output_v_path, mode="w", format='NETCDF4')
    output_v.createDimension('v_time', None)
    output_v.createDimension('v_lev', len(level))
    output_v.createDimension('v_lat', len(lat))
    output_v.createDimension('v_lon', len(lon))
    v_times= output_v.createVariable('v_time', 'i4', ('v_time',))
    v_levs = output_v.createVariable('v_lev', 'i4', ('v_lev',))
    v_lats = output_v.createVariable('v_lat', 'f4', ('v_lat',))
    v_lons = output_v.createVariable('v_lon', 'f4', ('v_lon',))
    v      = output_v.createVariable('v', 'f4', ('v_time','v_lev', 'v_lat', 'v_lon',))
    v_times[:]=np.array(time_v)
    v_levs[:] =np.array(level)
    v_lats[:] =np.array(lat)
    v_lons[:] =np.array(lon)
    vwind_time_lev=v_yeari.variables['v'][:,:,:,240:]
    v_time_lev_trends=time_i*v_trends
    v[:,:,:,:] = vwind_time_lev-v_time_lev_trends
    
    v_yeari.close()
    output_v.close()
    
#%%

q_trends=mk_results[1,:,:,:]
u_trends=mk_results[3,:,:,:]
v_trends=mk_results[5,:,:,:]

q_1950=nc.Dataset('I:/vairable_for_ar/humid/h1950.nc')
level = q_1950.variables['level']
lat   = q_1950.variables['latitude']
lon   = q_1950.variables['longitude'][240:]

for year_i in range(1950,2021,1):
    
    nc_minus_slope_q(year_i,q_trends,level,lat,lon)
    nc_minus_slope_u(year_i,u_trends,level,lat,lon)
    nc_minus_slope_v(year_i,v_trends,level,lat,lon)
    
    print(year_i)
    
#%% create IWV/IVT
# EXP0_1 remove q trends
# EXP0_2 remove winds trends
# EXP0_3 remove q and winds trends
def cal_ivt(i,humid,vwind,uwind,level):
    
    if exp=='EXP0_1':
        humid_time=humid.variables['q'][i,:,:,:]
        vwind_time=vwind.variables['v'][i,:,:,240:]
        uwind_time=uwind.variables['u'][i,:,:,240:]
    elif exp=='EXP0_2':
        humid_time=humid.variables['q'][i,:,:,240:]
        vwind_time=vwind.variables['v'][i,:,:,:]
        uwind_time=uwind.variables['u'][i,:,:,:]
    elif exp=='EXP0_3':
        humid_time=humid.variables['q'][i,:,:,:]
        vwind_time=vwind.variables['v'][i,:,:,:]
        uwind_time=uwind.variables['u'][i,:,:,:]
    h_v=humid_time*vwind_time
    h_u=humid_time*uwind_time
    ivt_v=100*simps(h_v,level,axis=0)/9.8
    ivt_u=100*simps(h_u,level,axis=0)/9.8
    ivt=np.sqrt(ivt_v*ivt_v+ivt_u*ivt_u)
    return ivt

def cal_iwv(i,humid,level):
    
    if exp=='EXP0_1':
        humid_time=humid.variables['q'][i,:,:,:]
    elif exp=='EXP0_2':
        humid_time=humid.variables['q'][i,:,:,240:]
    elif exp=='EXP0_3':
        humid_time=humid.variables['q'][i,:,:,:]
    iwv=100*simps(humid_time,level,axis=0)/9.8
    return iwv

#%% create IWV/IVT

for year_i in range(1950,2021,1):
    
    exp='EXP0_3'
    
    outputIVT_path = 'I:/AR_EXP/'+exp+'/IVT/IVT'+str(year_i)+'.nc';
    outputIWV_path = 'I:/AR_EXP/'+exp+'/IWV/IWV'+str(year_i)+'.nc';
    
    if exp=='EXP0_1':
        humid_path='G:/vairable_for_ar/humid_no_trends/h'+str(year_i)+'.nc';
        vwind_path='I:/vairable_for_ar/vwind/vwind'+str(year_i)+'.nc';
        uwind_path='I:/vairable_for_ar/uwind/uwind'+str(year_i)+'.nc';
        humid=nc.Dataset(humid_path)
        vwind=nc.Dataset(vwind_path)
        uwind=nc.Dataset(uwind_path)
        time=humid.variables['q_time']
        level=humid.variables['q_lev']
        lat=humid.variables['q_lat']
        lon=humid.variables['q_lon']
        
    elif exp=='EXP0_2':
        humid_path='I:/vairable_for_ar/humid/h'+str(year_i)+'.nc';
        vwind_path='G:/vairable_for_ar/vwind_no_trends/v'+str(year_i)+'.nc';
        uwind_path='G:/vairable_for_ar/uwind_no_trends/u'+str(year_i)+'.nc';
        humid=nc.Dataset(humid_path)
        vwind=nc.Dataset(vwind_path)
        uwind=nc.Dataset(uwind_path)
        time=humid.variables['time']
        level=humid.variables['level']
        lat=humid.variables['latitude']
        lon=humid.variables['longitude'][240:]
        
    elif exp=='EXP0_3':
        humid_path='G:/vairable_for_ar/humid_no_trends/h'+str(year_i)+'.nc';
        vwind_path='G:/vairable_for_ar/vwind_no_trends/v'+str(year_i)+'.nc';
        uwind_path='G:/vairable_for_ar/uwind_no_trends/u'+str(year_i)+'.nc';
        humid=nc.Dataset(humid_path)
        vwind=nc.Dataset(vwind_path)
        uwind=nc.Dataset(uwind_path)
        time=humid.variables['q_time']
        level=humid.variables['q_lev']
        lat=humid.variables['q_lat']
        lon=humid.variables['q_lon']

    outputIVT = nc.Dataset(outputIVT_path, mode="w", format='NETCDF4')
    Itime = outputIVT.createDimension('Itime', len(time))
    Ilat = outputIVT.createDimension('Ilat', len(lat))
    Ilon = outputIVT.createDimension('Ilon', len(lon))
    IVTtimes = outputIVT.createVariable('Itime', 'i4', ('Itime',))
    IVTlats = outputIVT.createVariable('Ilat', 'f8', ('Ilat',))
    IVTlons = outputIVT.createVariable('Ilon', 'f8', ('Ilon',))
    IVT_YEAR = outputIVT.createVariable('IVT_YEAR', 'f8', ('Itime', 'Ilat', 'Ilon',))
    IVTtimes[:]=np.array(time);
    IVTlats[:]=np.array(lat);
    IVTlons[:]=np.array(lon);
    
    outputIWV = nc.Dataset(outputIWV_path, mode="w", format='NETCDF4')
    IWVtime = outputIWV.createDimension('IWVtime', len(time))
    IWVlat = outputIWV.createDimension('IWVlat', len(lat))
    IWVlon = outputIWV.createDimension('IWVlon', len(lon))
    IWVtimes = outputIWV.createVariable('IWVtime', 'i4', ('IWVtime',))
    IWVlats = outputIWV.createVariable('IWVlat', 'f8', ('IWVlat',))
    IWVlons = outputIWV.createVariable('IWVlon', 'f8', ('IWVlon',))
    IWV_YEAR = outputIWV.createVariable('IWV_YEAR', 'f8', ('IWVtime', 'IWVlat', 'IWVlon',))
    IWVtimes[:]=np.array(time);
    IWVlats[:]=np.array(lat);
    IWVlons[:]=np.array(lon);
    
    for i in range(len(time)):
        
        IVT_YEAR[i,:,:]=cal_ivt(i,humid,vwind,uwind,level)
        IWV_YEAR[i,:,:]=cal_iwv(i,humid,level)
    
    humid.close()
    vwind.close()
    uwind.close()
    outputIVT.close()
    outputIWV.close()
    
    print(year_i)