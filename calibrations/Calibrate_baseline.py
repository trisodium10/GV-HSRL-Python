# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 12:47:00 2017

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
import datetime as datetime

import scipy
from scipy import optimize

import LidarProfileFunctions as lp
import GVHSRLlib as gv

"""
Write a baseline file with all zeros
"""

# input and raw_input are the same for python 3 and 2 respectively
# this makes it so input always accepts a string
try:
    input=raw_input
except NameError:
    pass

save_path_ez = os.path.abspath(__file__+'/../cal_files/')+'/'

#save_file_path = '/Users/mhayman/Documents/Python/Lidar/'
#save_file_path = '/h/eol/mhayman/HSRL/hsrl_processing/hsrl_configuration/projDir/calfiles/'
save_file_path = save_path_ez

basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path
#basepath = '/Users/mhayman/Documents/HSRL/GVHSRL_data/'  # local computer data path

tres = 4*60.0  # resolution in time in seconds (0.5 sec)

sg_win = 11
sg_order = 5
        
year_in = 2017
month_in = 10
day_in = 19
start_hr = 20
stop_hr = 2

print('This program only writes empty (zeros) baseline files')

print('Default Test Date:')
print('(M/D/Y) %d/%d/%d, starting %.1f UTC for %.1f h'%(month_in,day_in,year_in,start_hr,stop_hr))
if input('Run this default date? [y/n]') != 'y':
    print("Enter search range for diff_geo calibration:")
    year_in = np.int(input("Year: "))
    month_in = np.int(input("Month (#): "))
    day_in = np.int(input("Day: "))
    start_hr = np.float(input("Start Hour (UTC): "))
    stop_hr = np.float(input("Duration (hours): "))

cal_start = datetime.datetime(year_in,month_in,day_in)+datetime.timedelta(hours=start_hr)
cal_stop = cal_start + datetime.timedelta(hours=stop_hr)


var_1d_list = ['total_energy','RemoveLongI2Cell'\
    ,'TelescopeDirection','TelescopeLocked','polarization','DATA_shot_count']  # 'DATA_shot_count'


var_2d_list = ['cross','molecular','combined_hi','combined_lo']

[timeD,time_dt,time_sec],var_1d_data, profs = gv.load_raw_data(cal_start,cal_stop,var_2d_list,var_1d_list,basepath=basepath,verbose=True,as_prof=True)

master_time = np.arange(time_sec[0]-tres/2,time_sec[-1]+tres/2,tres)

t_profs,var_1d = gv.var_time_resample(master_time,time_sec,var_1d_data,average=True)

range_int = {}
plt.figure()
for ai, var in enumerate(profs.keys()):
    profs[var].time_resample(tedges=master_time,update=True,remainder=False)
    range_int[var] = np.nansum(profs[var].profile,axis=1)
    plt.plot(profs[var].time/3600.0,range_int[var])
    if ai == 0:
        cal_bool = np.logical_and(range_int[var] < 28000, range_int[var] > 2000)
    else:
        cal_bool = np.logical_and(cal_bool,np.logical_and(range_int[var] < 28000, range_int[var] > 2000))
plt.legend(list(range_int.keys()))
plt.grid(b=True)
plt.ylabel('Range Integrated Counts')
plt.xlabel('time [h-UTC]')
#    profs[var].bg_subtract(BGIndex)

cal_times = np.nonzero(cal_bool)[0]

#lp.plotprofiles(profs)
fig_data = lp.pcolor_profiles([profs['combined_hi'],profs['molecular']],ylimits=[0,12],climits=[[1e-4,1e4],[1e-4,1e4]]) 
fig_data[1][1].plot(profs['molecular'].time/3600,1+4*cal_bool,'r--')      
plt.show(block=False)

#t_profs = profs['molecular'].time.copy()
#if all(var_1d_data['RemoveLongI2Cell'] > 50):
#    print('No intervals found with I2 cell removed')
#    RunCal = False
#
#else:

RunCal = True
# find times where the i2 cell was inserted/removed
ical = np.nonzero(np.diff(1.0*cal_bool) != 0)[0]
if not cal_bool[0]:
    i0 = ical[::2]  # start indices
    i1 = ical[1::2] # stop indices
else:
    i0 = np.concatenate((np.zeros(1),ical[1::2]))  # start indices
    i1 = ical[::2] # stop indices
    
# ask user to select the interval to use
print('Found calibration intervals:')
for ai in range(len(i0)):
    if ai < i1.size:
        print('%d.)  %.2f - %.2f UTC (%.2f h)'%(ai,t_profs[i0[ai]]/3600,t_profs[i1[ai]]/3600,(t_profs[i1[ai]]-t_profs[i0[ai]])/3600))
    else:
        print('%d.)  %.2f - End of file UTC (%.2f h)'%(ai,t_profs[i0[ai]]/3600,(t_profs[-1]-t_profs[i0[ai]])/3600))
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
    i0 = np.array(list(i0) + np.argmin(np.abs(t_profs-t_input1)))
    i1 = np.array(list(i1) + np.argmin(np.abs(t_profs-t_input2)))
else:
    # the user selects a pre-defined calibration interval
    time_range = [1*60+t_profs[i0[cal_index]],t_profs[i1[cal_index]]-1*60]

if RunCal:
    # integrate the profiles over the selected time region
    for pname in profs.keys():
        profs[pname].slice_time(time_range)
        profs[pname].time_integrate(avg=True)

    avg_energy = np.nanmean(var_1d['total_energy'][i0[cal_index]:i1[cal_index]])
    lp.plotprofiles(profs)

    baseline_data = {}
    for var in profs.keys():
        # estimate baseline as a superposition of exponentials
        sol_list =[]
        error_list = []
        fit_sol_list = []
        mean_fit_error = []
        ifit= [3000,3999]
        
        fit_prof0 = profs[var].profile.flatten()
        fit_std0 = np.sqrt(profs[var].profile_variance.flatten())
        
        next_index = np.nonzero((fit_prof0[:ifit[1]-1]-fit_prof0[ifit[1]])/np.sqrt(fit_std0[ifit[1]]**2+fit_std0[:ifit[1]-1]**2)>4.0)[0]
        if len(next_index) > 0:
            ifit[0] = next_index[-1]
        
        fit_sol_sum = np.zeros(profs[var].profile.size)
        loop_on = True
        ai = 0
        while loop_on:
#            print('iteration: %d'%ai)
#            print('fit range: %d - %d'%tuple(ifit))
            
                
            fit_prof = fit_prof0[ifit[0]:ifit[1]]
            fit_range = profs[var].range_array[ifit[0]:ifit[1]]
            
            fit_prof = np.log(fit_prof)
            inan = np.nonzero(np.isnan(fit_prof))[0]
            fit_prof = np.delete(fit_prof,inan)
            fit_range = np.delete(fit_range,inan)
            
            mat = np.ones((fit_range.size,2))
            mat[:,1] = fit_range
            
            sol_list = sol_list+[np.dot(np.linalg.pinv(mat),fit_prof[:,np.newaxis])]
            error_list = error_list+[np.mean((fit_prof-np.dot(mat,sol_list[ai]))**2)]
            
            matsol = np.ones((profs[var].range_array.size,2))
            matsol[:,1] = profs[var].range_array
            fit_sol = np.dot(matsol,sol_list[ai])   
            
            fit_sol_sum = fit_sol_sum + np.exp(fit_sol.flatten())  
            mean_fit_error = mean_fit_error + [np.nanmean((profs[var].profile.flatten()[ifit[0]:]-fit_sol_sum[ifit[0]:])**2/fit_std0[ifit[0]:]**2)] 
            
            if ai > 0:
                if np.abs(mean_fit_error[ai]-mean_fit_error[ai-1])/mean_fit_error[ai-1] > np.max([1.0/ai,0.1]):
                    loop_on = False
            
            if loop_on:
                fit_sol_list = fit_sol_list+[fit_sol]
                fit_prof0 = fit_prof0-np.exp(fit_sol.flatten())
            
                ifit_last = ifit.copy()
           
                
                next_index = np.nonzero(np.abs(fit_prof0[:ifit[0]-1]/fit_std0[:ifit[0]-1]) > 4.0)[0]
                if len(next_index) > 0:
                    ifit[1] = next_index[-1]
                else:
                    loop_on = False
                
                next_index = np.nonzero(np.abs((fit_prof0[:ifit[1]-1]-fit_prof0[ifit[1]])/np.sqrt(fit_std0[ifit[1]]**2+fit_std0[:ifit[1]-1]**2))>4.0)[0]
                if len(next_index):
                    ifit[0] = next_index[-1]
                else:
                    loop_on = False
            
        ##        ifit_last = ifit.copy()
        #        ifit[1] = ifit[0]-1
        #        plt.figure();
        #        plt.plot(fit_prof0[:ifit[1]-1]/fit_std0[:ifit[1]-1])
        #        next_index = np.nonzero(np.abs(fit_prof0[:ifit[1]-1]/fit_std0[:ifit[1]-1]) > 2.0)[0]
        #        if len(next_index) > 0:
        #            ifit[0] = next_index[-1]
        #        else:
        #            loop_on = False
        #            
        #        if ifit[0] < 41:
        #            loop_on = False
                   
                if ai > 0:
                    if (mean_fit_error[ai]-mean_fit_error[ai-1])/mean_fit_error[ai-1] > np.max([1.0/ai,0.1]):
                        loop_on = False
            ai = ai + 1
                
                
                
            
        plt.figure()
        plt.plot(np.array(mean_fit_error))
        plt.title('fit error ' + var)
        plt.xlabel('function iteration')
        plt.ylabel('variance weighted error')
        plt.grid(b=True)
        
        fit_sol_sum = np.zeros(profs[var].profile.size)
        for ai in range(len(fit_sol_list)):
            fit_sol_sum = fit_sol_sum + np.exp(fit_sol_list[ai].flatten())

        imerge = np.nonzero(np.abs((fit_sol_sum[:ifit_last[1]]-profs[var].profile.flatten()[:ifit_last[1]])/fit_std0[:ifit_last[1]]) > 1.0)[0]+1
        sol_merge = np.concatenate((profs[var].profile.flatten()[:imerge[-1]],fit_sol_sum[imerge[-1]:]))
        
        plt.figure()
        plt.semilogy(profs[var].profile.flatten())
        plt.semilogy(fit_sol_sum)
        plt.semilogy(sol_merge,'r--')
        plt.legend(['raw signa','fit signal','merged'])
        plt.title('baseline '+var)
        plt.grid(b=True)
        
        plt.show(block=False)
        
        print('Select baseline data to use in file:')
        print('0.) Raw Data')
        print('1.) Merged Data')
        print('2.) All Zeros')
        data_sel = np.int(input(': '))
        if data_sel == 0:
            baseline_data[var] = profs[var].profile.flatten()
        elif data_sel == 2:
            baseline_data[var] = np.zeros(profs[var].profile.size)
        else:
            baseline_data[var] = sol_merge
          
    
    
save_cal = input("Save Calibrations[y/n]")
    
if save_cal == 'y' or save_cal == 'Y':

    save_file_name = 'baseline_correction_'+cal_start.strftime('%Y%m%dT%H%M')+'.blc'
    save_ez = 'baseline_'+cal_start.strftime('%Y%m%dT%H%M')+'_tmp'
    # write out zeros for baseline correction
    raw_data = np.zeros((4000,7))
    write_data = np.zeros((4000,7))
    if 'combined_hi' in profs.keys():
        write_data[:,1] = baseline_data['combined_hi']
        raw_data[:,1] = profs['combined_hi'].profile.flatten()
    if 'combined_lo' in profs.keys():
        write_data[:,2] = baseline_data['combined_lo']
        raw_data[:,2] = profs['combined_lo'].profile.flatten()
    if 'molecular' in profs.keys():
        write_data[:,3] = baseline_data['molecular']
        raw_data[:,3] = profs['molecular'].profile.flatten()
    if 'cross' in profs.keys():
        write_data[:,4] = baseline_data['cross']
        raw_data[:,4] = profs['cross'].profile.flatten()
    
    write_data[:,0] = np.arange(write_data.shape[0])
    raw_data[:,0] = np.arange(write_data.shape[0])
    
    header_str = 'profile data from ' + cal_start.strftime('%d-%b-%y %H:%M')+' -->'+cal_stop.strftime('%H:%M UTC') + '\n'
    header_str = header_str + 'ave energy per shot= %.5f    mJ\n'%(avg_energy*50e-9)
    header_str = header_str + 'bin_num\tcombined_hi\tcombined_lo\tmolecular\tcrosspol\tmol_I2A\t\tcomb_1064'
    
    np.savetxt(save_file_path+save_file_name,write_data,fmt='\t%d\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e',header=header_str,comments='#')
    
    np.savez(save_file_path+save_ez,write_data,raw_data = raw_data, \
        cal_range = np.array([time_dt[i0[cal_index]],time_dt[i1[cal_index]]]), \
        avg_energy = avg_energy)
