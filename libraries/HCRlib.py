#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 15:12:57 2019

@author: mhayman


Library for reading in HCR data


"""
import numpy as np
import matplotlib.pyplot as plt
import datetime


import LidarProfileFunctions as lp
import LidarPlotFunctions as lplt

import glob

import netCDF4 as nc4


hcr_data_defs = {'DBZ':
                    {'label':'dBZ',
                     'units':'dB',
                     'descript':'Log backscatter'},
                'LDR':
                    {'label':'LDR',
                     'units':'dB',
                     'descript':'Ratio of cross polarization to transmit polarization power'}
                    }
                    

def load_hcr_processed_files(time_start,time_stop,prof_hcr,data_path = None,verbose=True,date_reference=None,master_time=None):
    #time_start = datetime.datetime(2018,2,5,4,)
    #time_stop = datetime.datetime(2018,2,5,5,)
    
    #data_path = '/scr/rain1/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/'
    if data_path is None:
        #data_path = '/scr/snow2/rsfdata/projects/socrates/hcr/qc/cfradial/velcorr/2hz/'  # eventual permenant data path
        data_path = '/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/velcorr/2hz/'
    
    
#    prof_hcr = ['DBZ'] #,'DBZHC','DBZVC','VEL','WIDTH','SNR','LDR','LDRH','LDRV']
    hcr_t1d = ['time','azimuth','elevation','n_samples','latitude','longitude','altitude','ray_gate_spacing','roll','pitch']
#    hcr_r1d = ['range']
#    
#    verbose=True
#    master_time_edges = None
    
    hcr_tdata = {}
    hcr_profs = {}
    #hcr_remain_prof = {}
    if not master_time is None:
        regrid_time = True
        print('HCRlib load_hcr_processed_files()')
        print('   time grid matching is not implemented yet')
    
    
    dir_time = datetime.datetime(year=time_start.year,month=time_start.month,day=time_start.day)
    filefound = False # flag to indicate if a data file was ever loaded
    while dir_time <= time_stop:
        file_dir = data_path+dir_time.strftime('%Y%m%d')+'/'
        # typical file
        # cfrad.20180203_235900.119_to_20180204_000000.032_HCR_SUR.nc
        flist = sorted(glob.glob(file_dir+'cfrad.*_to_*_HCR_SOCRATES.nc'))
    
        for filename in flist:
            filestr = filename[len(file_dir):]  # remove preceeding directory to parse
            fsplit = filestr.split('.')
            file_start_time = datetime.datetime.strptime(fsplit[1],'%Y%m%d_%H%M%S')
            file_stop_time = datetime.datetime.strptime(fsplit[2].split('_to_')[1],'%Y%m%d_%H%M%S')
            
             # check if this file is in the search time period
            load_cond = (file_start_time >= time_start and file_start_time < time_stop) or \
                (file_stop_time < time_stop and file_stop_time > time_start) or \
                (file_start_time <= time_start) and (file_stop_time >= time_stop)
    
            if load_cond:
                filefound = True  # at least one file was found
                if verbose:
                    print( 'Loading '+file_start_time.strftime('%Y-%b-%d %H:%M:%S to ')+file_stop_time.strftime('%Y-%b-%d %H:%M:%S'))
                    print( '['+filestr+']')
                hcr_tdata=lp.load_nc_vars(filename,hcr_t1d,readfull=False,var_data0=hcr_tdata)  
                
    #            hcr_tdata0=lp.load_nc_vars(filename,['time','ray_gate_spacing'],readfull=False) # file grid data needed for profiles
                
                
                
                with nc4.Dataset(filename,'r') as f:
                    for pvar in prof_hcr:
                        loaded_prof = lp.LidarProfile(f.variables[pvar][:],f.variables['time'][:],\
                                label=hcr_data_defs[pvar]['label'],\
                                descript = hcr_data_defs[pvar]['descript'],\
                                lidar='HCR',StartDate=file_start_time,binwidth=np.median(f.variables['ray_gate_spacing'][:])*2/lp.c)
                        loaded_prof.range_array = f.variables['range'][:].copy()
                        loaded_prof.profile_type = hcr_data_defs[pvar]['units']
    #                    if regrid_time:
    #                        hcr_remain_profs[pvar] = loaded_prof.time_resample(tedges=master_time_edges,update=True,remainder=True,average=True)
                        
                        if pvar in hcr_profs.keys():
                            hcr_profs[pvar].cat_time(loaded_prof,front=False)
                        else:
                            hcr_profs[pvar] = loaded_prof.copy()
                            
        
        dir_time+=datetime.timedelta(days=1)
    if filefound:
        if not date_reference is None:
            for pvar in hcr_profs.keys():
                time_diff = hcr_profs[pvar].StartDate-date_reference
                hcr_profs[pvar].time += time_diff.total_seconds()
                hcr_profs[pvar].StartDate = date_reference
                
        hcr_tdata['time'] = hcr_profs[pvar].time.copy()  # change time variable so it doesn't reference the file start time
        
        hcr_tdata['radar_pointing'] = np.array([-np.cos(hcr_tdata['elevation']*np.pi/180)*np.cos(hcr_tdata['azimuth']*np.pi/180),
                                               -np.cos(hcr_tdata['elevation']*np.pi/180)*np.sin(hcr_tdata['azimuth']*np.pi/180),
                                               -np.sin(hcr_tdata['elevation']*np.pi/180)])
    else:
        # no file was found for requested time period
        hcr_profs = None
        hcr_tdata = None
        print('No HCR data found for '+time_start.strftime('%Y-%m-%d %H:%M:%S to ')+time_stop.strftime('%Y-%m-%d %H:%M:%S'))
        
    return hcr_profs,hcr_tdata