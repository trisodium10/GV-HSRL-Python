# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 14:50:44 2018

@author: mhayman
"""

import numpy as np
import matplotlib.pyplot as plt
import LidarProfileFunctions as lp
import LidarPlotFunctions as lplt
import datetime

import matplotlib.dates as mdates

#import os
#import sys
#
#import glob


import GVHSRLlib as gv

import copy


### Profile through aggrigates and plates RF03
#file_list = ['/h/eol/mhayman/HSRL/Plots/PTV/Denoise_20180123/'+'20180123_26_26_UTC_20180524_0851.nc']  # start at 12 minutes
## Profile through aggrigates and plates RF03
#file_list = ['/h/eol/mhayman/HSRL/Plots/PTV/Denoise_20180123/'+'20180123_26_26_UTC_20180604_1434.nc']  # start at 11 minutes

# Pristine columns in socrates RF03
#file_list = ['/h/eol/mhayman/HSRL/Plots/PTV/Denoise_20180123/'+'20180123_24_24_UTC_20180601_1158.nc']

# Oriented ice case with aircraft banking in RF07
#file_list = ['/h/eol/mhayman/HSRL/Plots/PTV/Denoise_20180131/'+'20180131_6_6_UTC_20180525_0952.nc']

path = '/h/eol/mhayman/HSRL/Plots/PTV/Denoise_20180123/'
file_list = [path+'20180123_26_26_UTC_20180607_1159.nc',
             path+'20180123_26_26_UTC_20180604_1434.nc',path+'20180123_26_26_UTC_20180608_1203.nc',
             path+'20180123_26_26_UTC_20180607_1031.nc',path+'20180123_26_26_UTC_20180607_1054.nc',
             path+'20180123_26_26_UTC_20180607_1109.nc',path+'20180123_26_26_UTC_20180607_1117.nc',
             path+'20180123_26_26_UTC_20180607_1142.nc',path+'20180123_26_26_UTC_20180607_1151.nc']

filename = 'SOCRATESrf03_GVHSRL_Denoised'

save_data = False
save_figs = False

#file_list  = ['/scr/sci/mhayman/SOCRATES/range_centered/test_data/SOCRATESrf03/SOCRATESrf03_GVHSRL_20180123T0000_20180123T0018.nc']
Project = 'SOCRATES '


prof_list = ['Denoised_Aerosol_Backscatter_Coefficient','Denoised_Particle_Depolarization','Denoised_Lidar_Ratio','Denoised_Aerosol_Extinction_Coefficient',
             'Aerosol_Backscatter_Coefficient']
             
             
plot_settings = {
    'Aerosol_Backscatter_Coefficient':{
        'climits':[1e-8,1e-3],
        'scale':'log',
	'colormap':'jet'
        },
    'Denoised_Aerosol_Extinction_Coefficient':{
        'climits':[1e-5,1e-2],
        'scale':'log',
	'colormap':'viridis'
        },
    'Denoised_Lidar_Ratio':{
        'climits':[0,50],
        'scale':'linear',
	'colormap':'jet'
        },
    'Denoised_Particle_Depolarization':{
        'climits':[0,0.7],
        'scale':'linear',
	'colormap':'viridis'
        },
    'Volume_Depolarization':{
        'climits':[0,1.0],
        'scale':'linear',
	'colormap':'jet'
        },
    'Denoised_Aerosol_Backscatter_Coefficient':{
        'climits':[1e-7,1e-3],
        'scale':'log',
	'colormap':'jet'
        },
    'Merged_Combined_Channel':{
        'climits':[1e-1,1e4],
        'scale':'log',
	'colormap':'jet'
        },
    }
    
profs,lidar_data,aircraft_data = gv.load_GVHSRL_processed_files([file_list[0]],prof_list,load_mask = True)

for ai in range(1,len(file_list)):
    profs_temp,lidar_data_temp,aircraft_data_temp = gv.load_GVHSRL_processed_files([file_list[ai]],prof_list,load_mask = True)
    
    
    
    # find the overlaping regions
    t_cat = np.concatenate((profs_temp['Denoised_Aerosol_Backscatter_Coefficient'].time,profs['Denoised_Aerosol_Backscatter_Coefficient'].time))
    tunique, tcounts = np.unique(t_cat, return_counts=True)
    #t_overlap_index = np.nonzero(tcounts==2)
    t_overlap = tunique[tcounts==2]
    t_start_overlap = t_overlap[0]
    t_end_overlap = t_overlap[-1]
    tindex_1 = np.argmin(np.abs(t_start_overlap-profs[prof_list[0]].time))
    tindex_2 = np.argmin(np.abs(t_end_overlap-profs_temp[prof_list[0]].time))
    
    lidar_data_new = {}
    for var in lidar_data.keys():
        if lidar_data[var].ndim == 1:
            lidar_data_new[var] = np.concatenate((lidar_data[var][:tindex_1],0.5*(lidar_data[var][tindex_1:]+lidar_data_temp[var][:tindex_2+1]),lidar_data_temp[var][tindex_2+1:]))
        elif lidar_data[var].ndim == 2:
            lidar_data_new[var] = np.hstack((lidar_data[var][:,:tindex_1],0.5*(lidar_data[var][:,tindex_1:]+lidar_data_temp[var][:,:tindex_2+1]),lidar_data_temp[var][:,tindex_2+1:]))
        
    lidar_data = copy.deepcopy(lidar_data_new)
    
    aircraft_data_new = {}
    for var in aircraft_data.keys():
        if aircraft_data[var].ndim == 1:
            aircraft_data_new[var] = np.concatenate((aircraft_data[var][:tindex_1],0.5*(aircraft_data[var][tindex_1:]+aircraft_data_temp[var][:tindex_2+1]),aircraft_data_temp[var][tindex_2+1:]))
        elif aircraft_data[var].ndim == 2:
            aircraft_data_new[var] = np.hstack((aircraft_data[var][:,:tindex_1],0.5*(aircraft_data[var][:,tindex_1:]+aircraft_data_temp[var][:,:tindex_2+1]),aircraft_data_temp[var][:,tindex_2+1:]))
         
    aircraft_data = copy.deepcopy(aircraft_data_new)
    
    # break each profile into the unique data and the overlapping region (4 profiles total)
    # average the overlapping regions the concatenate the three remaining profiles
    prof_new = {}
    for var in profs.keys():
        prof_new[var] = profs[var].copy()
        prof_new[var].slice_time_index(time_lim=[0,tindex_1])
        pm1 = profs[var].copy()
        pm1.slice_time_index(time_lim=[tindex_1+1,profs[var].time.size-2])
        pm2 = profs_temp[var].copy()
        pm2.slice_time_index(time_lim=[1,tindex_2-1])
        p2 = profs_temp[var].copy()
        p2.slice_time_index(time_lim=[tindex_2,profs_temp[var].time.size-1])
        
        pm1.profile = 0.5*(pm1.profile+pm2.profile)  
        
        prof_new[var].cat_time(pm1,front=False)
        prof_new[var].cat_time(p2,front=False)
        
        profs[var] = prof_new[var].copy()

profs['Denoised_Lidar_Ratio'].mask(profs['Denoised_Lidar_Ratio'].profile==20.0)

start_date = profs['Denoised_Aerosol_Backscatter_Coefficient'].StartDate+datetime.timedelta(microseconds=np.int(1e6*profs['Denoised_Aerosol_Backscatter_Coefficient'].time[0]))
stop_date = profs['Denoised_Aerosol_Backscatter_Coefficient'].StartDate+datetime.timedelta(microseconds=np.int(1e6*profs['Denoised_Aerosol_Backscatter_Coefficient'].time[-1]))

filedate = start_date.strftime('%Y%m%dT%H%M%S_')+stop_date.strftime('%Y%m%dT%H%M%S')

# var='Denoised_Aerosol_Backscatter_Coefficient'
for var in profs.keys(): 
    lplt.pcolor_airborne(profs[var],lidar_pointing = lidar_data['lidar_pointing'],lidar_alt=aircraft_data['GGALT'],
              climits=plot_settings[var]['climits'],scale=plot_settings[var]['scale'], #ylimits=[3.0,7.0],
              cmap=plot_settings[var]['colormap'],title_add ='SOCRATES RF03 ',plotAsDays=False,
              t_axis_scale=100.0,h_axis_scale=2.0,plot_date=True,
              minor_ticks=0,major_ticks=3.0/360)
    if save_figs:
        plt.savefig(path+filename+'_'+var+'_'+filedate+'.png',dpi=600)


if save_data:
    save_nc_name = path+filename+'_'+filedate+'.nc'
    
    for var in profs.keys():
        profs[var].write2nc(save_nc_name)
    for var in aircraft_data.keys():
        lp.write_var2nc(aircraft_data[var],var,save_nc_name)
    for var in lidar_data.keys():
        lp.write_var2nc(lidar_data[var],var,save_nc_name)

data_2D = np.load('/h/eol/mhayman/2DC/SOCRATESrf03_lower_regularizer__Denoised_94200_94800.npz')
data_cdp = np.load('/h/eol/mhayman/2DC/SOCRATESrf03h_Denoised_ACDP_94200_94799.npz')

var='Denoised_Lidar_Ratio'
x_time = mdates.date2num([datetime.datetime.fromordinal(profs[var].StartDate.toordinal()) \
                + datetime.timedelta(microseconds=np.int(sec)) for sec in profs[var].time*1e6])

x_time_2DC = mdates.date2num([datetime.datetime.fromordinal(profs[var].StartDate.toordinal()) \
                + datetime.timedelta(microseconds=np.int(sec)) for sec in data_2D['time'][data_2D['i0']:data_2D['i1']]*1e6])
                
x_time_cdp = mdates.date2num([datetime.datetime.fromordinal(profs[var].StartDate.toordinal()) \
                + datetime.timedelta(microseconds=np.int(sec)) for sec in data_cdp['time']*1e6])
                
cdp_counts = data_cdp['hist_denoise']
oap_counts =np.sum(data_2D['hist_denoise'],axis=1)
oap_rbar = np.sum(data_2D['hist_denoise']*data_2D['rpart'][np.newaxis,1:],axis=1)/oap_counts

#myFmt = mdates.DateFormatter('%H:%M:%S')
myFmt = mdates.DateFormatter('%M:%S')
ticks = 30  # tick locations in seconds
#ialt = [0,1,2,3]
ialt = [3,4,5]
plt.figure(figsize=(13,10))
plt.subplot(611); 
plt.semilogy(x_time,profs['Denoised_Aerosol_Backscatter_Coefficient'].profile[:,ialt],'b');
plt.gca().xaxis.set_major_formatter(myFmt); major_ticks = np.int(np.round(ticks))
plt.gca().xaxis.set_major_locator(mdates.SecondLocator(interval=major_ticks))
plt.gca().tick_params(axis='x', which='major', labelsize=8); 
plt.grid(b=True)
tset = plt.xlim()
plt.ylabel(r'$\beta_{a}$ $[m^{-1} sr^{-1}]$') 
plt.subplot(612); 
plt.plot(x_time,profs['Denoised_Particle_Depolarization'].profile[:,ialt],'g');
plt.gca().xaxis.set_major_formatter(myFmt); major_ticks = np.int(np.round(ticks))
plt.gca().xaxis.set_major_locator(mdates.SecondLocator(interval=major_ticks))
plt.gca().tick_params(axis='x', which='major', labelsize=8); 
plt.grid(b=True)
plt.ylabel('$d_{a}$') 
plt.xlim(tset)
plt.subplot(613); 
plt.plot(x_time,profs['Denoised_Lidar_Ratio'].profile[:,ialt],'r')
plt.gca().xaxis.set_major_formatter(myFmt); major_ticks = np.int(np.round(ticks))
plt.gca().xaxis.set_major_locator(mdates.SecondLocator(interval=major_ticks))
plt.gca().tick_params(axis='x', which='major', labelsize=8); 
plt.grid(b=True)
plt.ylabel('$s_{LR}$ $[sr]$')
plt.xlim(tset)
plt.subplot(614); 
plt.plot(x_time_cdp,cdp_counts,'k')
plt.gca().xaxis.set_major_formatter(myFmt); major_ticks = np.int(np.round(ticks))
plt.gca().xaxis.set_major_locator(mdates.SecondLocator(interval=major_ticks))
plt.gca().tick_params(axis='x', which='major', labelsize=8); 
plt.grid(b=True)
plt.ylabel('CDP Counts')
plt.xlim(tset)
plt.subplot(615); 
plt.plot(x_time_2DC,oap_rbar,'m')
plt.gca().xaxis.set_major_formatter(myFmt); major_ticks = np.int(np.round(ticks))
plt.gca().xaxis.set_major_locator(mdates.SecondLocator(interval=major_ticks))
plt.gca().tick_params(axis='x', which='major', labelsize=8); 
plt.grid(b=True)
plt.ylabel('2DC Mean Radius')
plt.xlim(tset)
plt.subplot(616); 
plt.semilogy(x_time_2DC,oap_counts,'c')
plt.gca().xaxis.set_major_formatter(myFmt); major_ticks = np.int(np.round(ticks))
plt.gca().xaxis.set_major_locator(mdates.SecondLocator(interval=major_ticks))
plt.gca().tick_params(axis='x', which='major', labelsize=8); 
plt.grid(b=True)
plt.ylabel('2DC Counts')
plt.xlim(tset)
