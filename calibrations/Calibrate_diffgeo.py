# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 09:25:01 2017

@author: mhayman
"""


import sys
import os

# add the path to GVHSRLlib manually
library_path = os.path.abspath(__file__+'/../../libraries/')
print(library_path)
if library_path not in sys.path:
    sys.path.append(library_path)


import numpy as np
import matplotlib.pyplot as plt
#from scipy.io import netcdf
#import netCDF4 as nc4
import LidarProfileFunctions as lp
#import scipy.interpolate
import datetime
#import glob

import GVHSRLlib as gv
    


# input and raw_input are the same for python 3 and 2 respectively
# this makes it so input always accepts a string
try:
    input=raw_input
except NameError:
    pass



save_path_ez = os.path.abspath(__file__+'/../cal_files/')+'/'

save_file_path = save_path_ez
#save_file_path = '/Users/mhayman/Documents/Python/Lidar/'
#save_file_path = '/h/eol/mhayman/HSRL/hsrl_processing/hsrl_configuration/projDir/calfiles/'

settings = {
    'deadtime_correct':True  # correct deadtime
    }

airborne = False

sg_win = 11
sg_order = 5
        
year_in = 2017
month_in = 10
day_in = 25
start_hr = 18
stop_hr = 4

print('Default Test Date:')
print('(M/D/Y) %d/%d/%d, starting %.1f UTC for %.1f h'%(month_in,day_in,year_in,start_hr,stop_hr))
if input('Run this default date? [y/n]') != 'y':
    print("Enter search range for diff_geo calibration:")
    year_in = np.int(input("Year: "))
    month_in = np.int(input("Month (#): "))
    day_in = np.int(input("Day: "))
    start_hr = np.float(input("Start Hour (UTC): "))
    stop_hr = np.float(input("Duration (hours): "))
    airborne_str = input("Is the lidar airborne? [y/n]")
    
    if airborne_str == 'y' or airborne_str == 'Y':
        airborne = True

cal_start = datetime.datetime(year_in,month_in,day_in)+datetime.timedelta(hours=start_hr)
cal_stop = cal_start + datetime.timedelta(hours=stop_hr)



#Day = 6
#Month = 10
#Year = 2017
#HourLim = np.array([0,24])  # Limits on the processing time

tres = 1*60.0  # resolution in time in seconds (0.5 sec)
zres = 10.0  # resolution in altitude points (7.5 m)

# index for where to treat the profile as background only
BGIndex = -100; # negative number provides an index from the end of the array
platform = 'ground' # 'ground' or 'airborne'.  If 'airborne' it needs an aircraft netcdf.
MaxAlt = 10e3


pol_xtalk = 0.015

#kB = 1.3806504e-23;
#c = 3e8


FilterI2 = False  # only include data where the I2 cell is removed

# list of 1D variables to load
var_1d_list = ['total_energy','RemoveLongI2Cell'\
    ,'TelescopeDirection','TelescopeLocked','polarization','DATA_shot_count']  # 'DATA_shot_count'

# list of 2D variables (profiles) to load
var_2d_list = ['molecular','combined_hi','combined_lo']  # 'cross',

import json
filename = __file__
cal_file_path = os.path.abspath(filename+'/../cal_files/')+'/'
cal_file = cal_file_path + 'gv_calvals.json'

with open(cal_file,"r") as f:
    cal_json = json.loads(f.read())
f.close()

dead_time_list = lp.get_calval(cal_start,cal_json,"Dead_Time",returnlist=var_2d_list)
dead_time = dict(zip(var_2d_list,dead_time_list))


#basepath = '/scr/eldora1/HSRL_data/'  # old path - still works with link from HSRL_data to /hsrl/raw/
#basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path
basepath = '/scr/rain1/rsfdata/projects/socrates/hsrl/raw/'  # SOCRATES Path
#basepath = '/Users/mhayman/Documents/HSRL/GVHSRL_data/'  # local computer data path


# grab raw data from netcdf files
[timeD,time_dt,time_sec],var_1d_data, profs = gv.load_raw_data(cal_start,cal_stop,var_2d_list,var_1d_list,basepath=basepath,verbose=True,as_prof=True)

# plot I2 cell status
#plt.figure(); 
#plt.plot(time_sec/3600,var_1d_data['RemoveLongI2Cell'])
#plt.grid(b=True)
#plt.xlabel('time [h-UTC]')
#plt.ylabel('Value')
#plt.title('RemoveLongI2Cell')

# equation for estimating transmit power is assumed based on known approximate
# DAQ output.  The conversion used here could be incorrect.
fig, ax1 = plt.subplots(); 
ax1.plot(time_sec/3600.0,0.5*var_1d_data['total_energy']/var_1d_data['DATA_shot_count'][:,0],'b'); 
ax1.set_xlabel('time [h-UTC]'); 
ax1.set_ylabel('Transmit Power [mW]', color='b')
ax1.tick_params('y', colors='b')
ax1.grid(b=True)
ax2 = ax1.twinx(); 
ax2.plot(time_sec/3600,var_1d_data['RemoveLongI2Cell'],'r');
ax2.set_ylabel('RemoveLongI2Cell', color='r')
ax2.tick_params('y', colors='r')


# set the master time to match all 2D profiles to
# (1d data will not be resampled)
master_time = np.arange(time_sec[0]-tres/2,time_sec[-1]+tres/2,tres)

for var in profs.keys():
    if settings['deadtime_correct']:
        if hasattr(profs[var],'NumProfsList'):
            profs[var].nonlinear_correct(dead_time[var],laser_shot_count=2000*profs[var].NumProfsList[:,np.newaxis],std_deadtime=5e-9)
        else:
            # number of laser shots is based on an assumption that there is one 0.5 second profile per time bin
            profs[var].nonlinear_correct(dead_time[var],laser_shot_count=2000,std_deadtime=5e-9)
    profs[var].time_resample(tedges=master_time,update=True,remainder=False)
    profs[var].bg_subtract(BGIndex)

## option to remove all data where the I2 cell is not removed.
## depricated
#if FilterI2:
#    i2_size = var_1d_data['RemoveLongI2Cell'].size  # size of array.  Don't apply to any arrays that don't have matching time dimensions
#    i2_rem = np.nonzero(var_1d_data['RemoveLongI2Cell'] < 50)[0]  # data points where the I2 cell is removed
#    
#    for var in var_1d_data.keys():
#            if var_1d_data[var].size == i2_size:
#                var_1d_data[var] = var_1d_data[var][i2_rem]
#            else:
#                print(var+' does not have time dimension that matches RemoveLongI2Cell')
#    
#    time_dt = time_dt[i2_rem]
#    time_sec = time_sec[i2_rem]

lp.plotprofiles(profs)
# overlay profiles and I2 cell status to help determine the calibration interval
fig_data = lp.pcolor_profiles([profs['combined_hi'],profs['molecular']],ylimits=[0,12],climits=[[1e-4,1e4],[1e-4,1e4]])  
fig_data[1][1].plot(time_sec/3600,var_1d_data['RemoveLongI2Cell']/25,'b--')      
plt.show(block=False)

# check if the I2 cell is ever removed in the data.  If not, kill the 
# process.  No point in running the calibration
if all(var_1d_data['RemoveLongI2Cell'] > 50):
    print('No intervals found with I2 cell removed')
    RunCal = False

else:
    RunCal = True
    # find times where the i2 cell was inserted/removed
    ical = np.nonzero(np.diff(var_1d_data['RemoveLongI2Cell'])!=0)[0]
    if var_1d_data['RemoveLongI2Cell'][0] > 50:
        i0 = ical[::2]  # start indices
        i1 = ical[1::2] # stop indices
    else:
        i0 = np.concatenate((np.zeros(1),ical[1::2]))  # start indices
        i1 = ical[::2] # stop indices
        
    # ask user to select the interval to use
    print('Found calibration intervals:')
    for ai in range(len(i0)):
        if ai < i1.size:
            print('%d.)  %.2f - %.2f UTC'%(ai,time_sec[i0[ai]]/3600,time_sec[i1[ai]]/3600))
        else:
            print('%d.)  %.2f - End of file UTC'%(ai,time_sec[i0[ai]]/3600))
            i1 = np.concatenate((i1,np.array([-1])))
    
    print('%d.)  custom'%(ai+1))    
    cal_index = np.int(input('Select Interval (invalid number quits cal): '))
    
        
    
    if cal_index > len(i0) or cal_index < 0:
        # if the user gives an out of bounds number, quit the routine
        RunCal = False
    elif cal_index == len(i0):
        # the user wants to enter a custom range
        t_input1 = np.float(input('Enter start time in h-UTC: '))
        t_input2 = np.float(input('Enter stop time in h-UTC:'))
        time_range = [t_input1, t_input2]
        i0 = np.array(list(i0) + np.argmin(np.abs(time_sec-t_input1)))
        i1 = np.array(list(i1) + np.argmin(np.abs(time_sec-t_input2)))
    else:
        # the user selects a pre-defined calibration interval
        time_range = [5*60+time_sec[i0[cal_index]],time_sec[i1[cal_index]]-5*60]
if RunCal:
    # integrate the profiles over the selected time region
    for pname in profs.keys():
        profs[pname].slice_time(time_range)
        profs[pname].time_integrate()
    
    # ratio the profiles to obtain the diff-geo corrections
    pHi = profs['combined_hi'].copy()   
    pLo = profs['combined_lo'].copy() 
#    pHiLo = profs['combined_hi'].copy()
    
    # using LidarProfile data types helps propagate the error
    pHi.divide_prof(profs['molecular'])   
    pLo.divide_prof(profs['molecular'])  
#    pHiLo.divide_prof(profs['combined_lo'])  
  
    hi_smooth = gv.savitzky_golay(pHi.profile.flatten(), sg_win, sg_order, deriv=0)  
    #    dhi_smooth = savitzky_golay(pHi.profile.flatten(), 11, 4, deriv=1)
        
    plt.figure()
    plt.plot(pHi.profile.flatten())
    plt.plot(hi_smooth)
    plt.title('Hi/Mol Ratio')
    plt.show(block=False)
    

    i_norm = 1000

    i_const = np.int(input('Make constant above index (e.g. 500): '))
    i_const_max = np.nonzero(pHi.SNR().flatten() < 0.2*pHi.SNR().flatten()[i_const])[0]
    i_const_max = i_const_max[np.nonzero(i_const_max > i_const)[0][0]]
    
    
    hi_diff_geo = hi_smooth
#    hi_diff_geo[i_const:] = hi_diff_geo[i_const]
#    plt.plot(hi_diff_geo)
    
    if airborne:
        hi_diff_const = np.nanmean(hi_diff_geo[i_const-4:i_const])
        hi_diff_geo[i_const:] = hi_diff_const
    else:
        hi_diff_geo = gv.fit_high_range(hi_diff_geo,pHi.profile_variance,i_const,i_const_max)
    plt.plot(hi_diff_geo,'--')
    
    hi_norm = hi_diff_geo[i_norm]
    hi_diff_geo = hi_diff_geo/hi_norm
    
    lo_smooth = gv.savitzky_golay(pLo.profile.flatten(), sg_win, sg_order, deriv=0) 
    
    plt.figure()
    plt.plot(pLo.profile.flatten())
    plt.plot(lo_smooth)
    plt.title('Lo/Mol Ratio')
    plt.show(block=False)
    
    i_const = np.int(input('Make constant above index (e.g. 590): '))
    i_const_max = np.nonzero(pLo.SNR().flatten() < 0.2*pLo.SNR().flatten()[i_const])[0]
    i_const_max = i_const_max[np.nonzero(i_const_max > i_const)[0][0]]
    
    
    lo_diff_geo = lo_smooth
#    lo_diff_geo[i_const:] = lo_diff_geo[i_const]
#    plt.plot(lo_diff_geo)
    
    if airborne:
        lo_diff_const = np.nanmean(lo_diff_geo[i_const-4:i_const])
        lo_diff_geo[i_const:] = lo_diff_const
    else:
        lo_diff_geo = gv.fit_high_range(lo_diff_geo,pLo.profile_variance,i_const,i_const_max)
    plt.plot(lo_diff_geo,'--')
    plt.show(block=False)
    
    lo_norm = lo_diff_geo[i_norm]
    lo_diff_geo = lo_diff_geo/lo_norm
    
    save_cal = input("Save Calibrations[y/n]")
    
    if save_cal == 'y' or save_cal == 'Y':
        write_data = np.ones((hi_diff_geo.size,4))
        write_data[:,0] = np.arange(hi_diff_geo.size)
        write_data[:,1] = hi_diff_geo
        write_data[:,2] = lo_diff_geo
        
        header_str = 'profile data from ' + time_dt[i0[cal_index]].strftime('%d-%b-%y %H:%M')+' -->'+time_dt[i1[cal_index]].strftime('%H:%M UTC') +'\n'
        header_str = header_str+'Gains normalized range bin = %d'%i_norm + '\n'
        header_str = header_str+'bin #\tdiff_hi/mol\tdiff_lo/mol\tdiff_i2a/mol'
    
        save_cal_file = 'diff_default_geofile_'+time_dt[i0[cal_index]].strftime('%Y%m%dT%H%M')+'.geo'

        save_cal_file_ez = 'diff_geo_GVHSRL'+time_dt[i0[cal_index]].strftime('%Y%m%d')
        
        print('Saving to:\n  '+save_file_path+save_cal_file)        
        
        np.savetxt(save_file_path+save_cal_file,write_data,fmt='\t%d\t%.4e\t%.4e\t%.4e',header=header_str)
        
        print('Saving python vars to:\n  '+save_path_ez+save_cal_file_ez+'.npz') 
        
        np.savez(save_path_ez+save_cal_file_ez,start_time=time_dt[i0[cal_index]],stop_time=time_dt[i1[cal_index]], \
            hi_diff_geo=hi_diff_geo,lo_diff_geo=lo_diff_geo, \
            hi_diff_var=pHi.profile_variance.data.flatten(),lo_diff_var=pLo.profile_variance.data.flatten(),\
            hi_prof = pHi.profile.data.flatten(),lo_prof = pLo.profile.data.flatten(),\
            sg_win = sg_win, sg_order = sg_order, i_norm = i_norm, 
            range_array=pHi.range_array,lo_norm=lo_norm,hi_norm=hi_norm, \
            TelescopeDirection = np.nanmean(var_1d_data['TelescopeDirection'][i0[cal_index]:i1[cal_index]]), \
            airborne=airborne)
        # note that if the lidar is airborne, it is probably a downward pointing calibration
        

