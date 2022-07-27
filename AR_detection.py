# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 20:06:47 2021

@author: Qiang Wang
"""

import pandas  as pd
import netCDF4 as nc
import numpy   as np
import sympy   as sp
from geopy.distance  import distance
from skimage.measure import label, regionprops
from scipy.io        import loadmat

#%% AR detection using IVT and IWV for EXP 1-3

def AR_detection_for_exp(exp):
    
    map_path1='I:/AR_EXP/'+exp+'_IVT_85TH.mat'
    map_path2='I:/AR_EXP/'+exp+'_IWV_85TH.mat'
    ivt_mat_85 = loadmat(map_path1)         #revise
    iwv_mat_85 = loadmat(map_path2)         #revise
    ivt_85th = np.array(pd.DataFrame(ivt_mat_85["var_85TH"])).T  
    iwv_85th = np.array(pd.DataFrame(iwv_mat_85["var_85TH"])).T  
    
    for Time_year in range(1950,2021,1):
    
        ivt_path='I:/AR_EXP/'+exp+'/IVT/IVT'+str(Time_year)+'.nc'
        iwv_path='I:/AR_EXP/'+exp+'/IWV/IWV'+str(Time_year)+'.nc'
        outputAR_path='I:/AR_EXP/'+exp+'/AR_records/AR_'+str(Time_year)+'.nc'
        
        ivt_year=nc.Dataset(ivt_path);
        iwv_year=nc.Dataset(iwv_path);
    
        time=ivt_year.variables['Itime'];
        lat =ivt_year.variables['Ilat'];
        lon =ivt_year.variables['Ilon'];
        ivt_all=ivt_year.variables['IVT_YEAR']
        iwv_all=iwv_year.variables['IWV_YEAR']
    
        outputAR = nc.Dataset(outputAR_path, mode="w", format='NETCDF4')
        outputAR.createDimension('ARtime', len(time))
        outputAR.createDimension('ARlat', len(lat))
        outputAR.createDimension('ARlon', len(lon))
        times = outputAR.createVariable('ARtime', 'i4', ('ARtime',))
        lats = outputAR.createVariable('ARlat', 'f8', ('ARlat',))
        lons = outputAR.createVariable('ARlon', 'f8', ('ARlon',))
        ARrecords = outputAR.createVariable('ARrecords', 'i1', ('ARtime', 'ARlat', 'ARlon',))
        times[:]=np.array(time);
        lats[:]=np.array(lat);
        lons[:]=np.array(lon);
        
        for i in range(len(time)):
    
            ivt=ivt_all[i,:,:]
            iwv=iwv_all[i,:,:]
            
            ivt[ivt<250]=0
            iwv[iwv<15]=0
            
            ivt_diff=ivt-ivt_85th
            ivt_diff[ivt_diff<=0]=0
            ivt_diff[ivt_diff>0]=1
            iwv_diff=iwv-iwv_85th
            iwv_diff[iwv_diff<=0]=0
            iwv_diff[iwv_diff>0]=1
            
            ar_timei= ivt_diff*iwv_diff
            
            # l>1500, r>2
            label_ar = label(ar_timei);
            ar_all=regionprops(label_ar);
            
            for regioni in ar_all:
                
                ari_centriod=regioni.centroid
                ari_grid=regioni.coords
                ari_diameter=regioni.major_axis_length
                ari_orien=regioni.orientation
                ari_eccentricity=regioni.eccentricity

                a1=ari_centriod[0]
                b1=ari_centriod[1]
                a2=a1 + 0.5*ari_diameter*sp.cos(ari_orien)
                b2=b1 + 0.5*ari_diameter*sp.sin(ari_orien)
                
                if round(b2)>=len(lon):
                    b2=b2-len(lon);
                if round(a2)>=len(lat):
                    a2=a2-len(lat)
                    
                lat1=lat[round(a1)]
                lon1=lon[round(b1)]
                lat2=lat[round(a2)]
                lon2=lon[round(b2)]
                
                ari_lenth = 2 * distance((lat1,lon1),(lat2,lon2)).km
                
                if ari_lenth<1500:
                    ar_timei[ari_grid[:,0],ari_grid[:,1]]=0
                if ari_eccentricity<0.866:
                    ar_timei[ari_grid[:,0],ari_grid[:,1]]=0
                    
            ARrecords[i,:,:]=ar_timei
        
        ivt_year.close()
        iwv_year.close()
        outputAR.close()
        
        print(exp,Time_year)

#%%ls
AR_detection_for_exp('XXX')
