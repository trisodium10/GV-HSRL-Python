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

#import scipy.interpolate
import datetime
#import glob

import GVHSRLlib as gv
import ExternalDataFunctions as ex
import LidarProfileFunctions as lp

import json

#import scipy as sp
    


# input and raw_input are the same for python 3 and 2 respectively
# this makes it so input always accepts a string
try:
    input=raw_input
except NameError:
    pass

#diff_geo_correct = True
#hsrl_rb_adjust = True
#RemoveCals = True # removes instances where the I2 cell is removed
#load_reanalysis = True
#use_bs_thr = True

cal_file_path = os.path.abspath(__file__+'/../../calibrations/cal_files/')+'/'
cal_file = cal_file_path + 'gv_calvals.json'

reanalysis_path = os.path.abspath(__file__+'/../../../external_data/')+'/'

save_path_ez = os.path.abspath(__file__+'/../cal_files/')+'/'

save_file_path = save_path_ez
#save_file_path = '/Users/mhayman/Documents/Python/Lidar/'
#save_file_path = '/h/eol/mhayman/HSRL/hsrl_processing/hsrl_configuration/projDir/calfiles/'

bs_th = 3e-6  # BS coeff threshold where profile isn't included in the geo overlap analysis

sg_win = 11
sg_order = 5
        
year_in = 2017
month_in = 10
day_in = 25
start_hr = 22
stop_hr = 2

print('Default Test Date:')
print('(M/D/Y) %d/%d/%d, starting %.1f UTC for %.1f h'%(month_in,day_in,year_in,start_hr,stop_hr))
if input('Run this default date? [y/n]') != 'y':
    print("Enter search range for geo calibration:")
    year_in = np.int(input("Year: "))
    month_in = np.int(input("Month (#): "))
    day_in = np.int(input("Day: "))
    start_hr = np.float(input("Start Hour (UTC): "))
    stop_hr = np.float(input("Duration (hours): "))

cal_start = datetime.datetime(year_in,month_in,day_in)+datetime.timedelta(hours=start_hr)
cal_stop = cal_start + datetime.timedelta(hours=stop_hr)



#Day = 6
#Month = 10
#Year = 2017
#HourLim = np.array([0,24])  # Limits on the processing time

settings = {
    'tres':1*60.0,  # resolution in time in seconds (0.5 sec)
    'zres':10.0,  # resolution in altitude points (7.5 m)
    
    # index for where to treat the profile as background only
    'BGIndex':-100, # negative number provides an index from the end of the array
    'platform':'ground', # 'ground' or 'airborne'.  If 'airborne' it needs an aircraft netcdf.
    'MaxAlt':40e3,
    
    'LRassumed':30,  # assumed lidar ratio for OD estimate
    
    
    'pol_xtalk':0.015,
    
    
    'FilterI2':False,  # only include data where the I2 cell is removed
    
    'diff_geo_correct':True,
    'hsrl_rb_adjust':True,
    'RemoveCals':True, # removes instances where the I2 cell is removed
    'load_reanalysis':True,
    'use_bs_thr':True,
    
    'EstimateExtinction':True,  # attempt to apply an extinction correction
    'ExtinctionAlt':5e3     # max altitude to apply extinction correction

    
    }

## list of 1D variables to load
#var_1d_list = ['total_energy','RemoveLongI2Cell'\
#    ,'TelescopeDirection','TelescopeLocked','polarization','DATA_shot_count']  # 'DATA_shot_count'
#
## list of 2D variables (profiles) to load
#var_2d_list = ['molecular','combined_hi','combined_lo']  # 'cross',

# list of 1D variables to load
var_1d_list = ['total_energy','RemoveLongI2Cell'\
    ,'TelescopeDirection','TelescopeLocked','polarization','DATA_shot_count']  # 'DATA_shot_count'

# list of 2D variables (profiles) to load
var_2d_list = ['molecular','combined_hi','combined_lo','cross']

#basepath = '/scr/eldora1/HSRL_data/'  # old path - still works with link from HSRL_data to /hsrl/raw/
#basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path
basepath = '/Users/mhayman/Documents/HSRL/GVHSRL_data/'  # local computer data path


# grab raw data from netcdf files
[timeD,time_dt,time_sec],var_1d_data, profs = gv.load_raw_data(cal_start,cal_stop,var_2d_list,var_1d_list,basepath=basepath,verbose=True,as_prof=True)
time_start = datetime.datetime(year_in,month_in,day_in)+datetime.timedelta(hours=start_hr)
time_stop = time_start + datetime.timedelta(hours=stop_hr)


## grab raw data from netcdf files
#[timeD,time_dt,time_sec],var_1d_data, profs = gv.load_raw_data(time_start,time_stop,var_2d_list,var_1d_list,basepath=basepath,verbose=True,as_prof=True)

# find instances in raw data where I2 cell is removed
if 'RemoveLongI2Cell' in var_1d_data.keys():
    cal_indices = np.nonzero(var_1d_data['RemoveLongI2Cell'] < 50)[0]
else:
    cal_indices = []

with open(cal_file,"r") as f:
    cal_json = json.loads(f.read())
f.close()
if settings['hsrl_rb_adjust']:
    mol_gain,diff_geo_file = lp.get_calval(time_start,cal_json,'Molecular Gain',cond=[['RB_Corrected','=','True']],returnlist=['value','diff_geo'])  
else:
    mol_gain,diff_geo_file = lp.get_calval(time_start,cal_json,"Molecular Gain",returnlist=['value','diff_geo'])
baseline_file = lp.get_calval(time_start,cal_json,"Baseline File")[0]
i2_file = lp.get_calval(time_start,cal_json,"I2 Scan")[0]
# load differential overlap correction
diff_data = np.load(cal_file_path+diff_geo_file)
baseline_data = np.load(cal_file_path+baseline_file)
i2_data = np.load(cal_file_path+i2_file)
#diff_data = np.load(cal_file_path+'diff_geo_GVHSRL20171025_tmp.npz')
#baseline_data = np.load(cal_file_path+'diff_geo_GVHSRL20171025_tmp.npz')

lidar_location = lp.get_calval(time_start,cal_json,"Location",returnlist=['latitude','longitude'])

tres = settings['tres']
BGIndex = settings['BGIndex']
MaxAlt = settings['MaxAlt']
# set the master time to match all 2D profiles to
# (1d data will not be resampled)
master_time = np.arange(time_sec[0]-tres/2,time_sec[-1]+tres/2,tres)

time_1d,var_1d = gv.var_time_resample(master_time,time_sec,var_1d_data,average=True)
int_profs = {}  # obtain time integrated profiles
for var in profs.keys():
    if settings['RemoveCals']:
        # remove instances where the I2 cell is removed
        profs[var].remove_time_indices(cal_indices)
    profs[var].time_resample(tedges=master_time,update=True,remainder=False)
#    int_profs[var] = profs[var].copy()
#    int_profs[var].time_integrate()
    
    if var == 'molecular':
        MolRaw = profs['molecular'].copy()
    
    profs[var].bg_subtract(BGIndex)
    if var == 'combined_hi' and settings['diff_geo_correct']:
        CombRaw = profs['molecular'].copy()
        profs[var].diff_geo_overlap_correct(diff_data['hi_diff_geo'])
    elif var == 'combined_lo' and settings['diff_geo_correct']:
        profs[var].diff_geo_overlap_correct(diff_data['lo_diff_geo'])
        profs[var].gain_scale(1.0/diff_data['lo_norm'])
    
#    profs[var].slice_range(range_lim=[0,MaxAlt])
#    int_profs[var].bg_subtract(BGIndex)
#    int_profs[var].slice_range(range_lim=[0,MaxAlt])




if settings['load_reanalysis']:
    pres,temp = ex.load_fixed_point_NCEP_TandP(profs['molecular'],lidar_location,reanalysis_path)
    beta_m = lp.get_beta_m(temp,pres,profs['molecular'].wavelength)
    
    
if settings['hsrl_rb_adjust']:
    print('Obtaining Rayleigh-Brillouin Correction')
    dnu = 20e6  # resolution
    nu_max = 10e9 # max frequency relative to line center
    nu = np.arange(-nu_max,nu_max,dnu)
    Ti2 = np.interp(nu,i2_data['freq']*1e9,i2_data['mol_scan'])  # molecular transmission
    
    Tc2 = np.interp(nu,i2_data['freq']*1e9,i2_data['combined_scan'])  # combined transmission
#    Trb = spec.RubidiumCellTransmission(nu+lp.c/Molecular.wavelength,RbCellTemp,RbCellPressure,RbCellLength,iso87=RbCellPurity)
#    Tm = SurfaceTemp_HSRL.profile-0.0065*Molecular.range_array[np.newaxis,:]
#    Pm = SurfacePres_HSRL.profile*(SurfaceTemp_HSRL.profile/Tm)**(-5.5)# /9.86923e-6  # give pressure in Pa
    
    beta_mol_norm = lp.RB_Spectrum(temp.profile.flatten(),pres.profile.flatten()*9.86923e-6,profs['molecular'].wavelength,nu=nu,norm=True)
    eta_i2 = np.sum(Ti2[:,np.newaxis]*beta_mol_norm,axis=0)
    eta_i2 = eta_i2.reshape(temp.profile.shape)
    profs['molecular'].multiply_piecewise(1.0/eta_i2)
    profs['molecular'].gain_scale(mol_gain)
    
    eta_c = np.sum(Tc2[:,np.newaxis]*beta_mol_norm,axis=0)
    eta_c = eta_c.reshape(temp.profile.shape)
    profs['combined_hi'].multiply_piecewise(1.0/eta_c)
    
#    if Denoise_Mol:
#        MolDenoise.multiply_piecewise(1.0/eta_i2)
#        MolDenoise.gain_scale(mol_gain)
else:
    # Rescale molecular channel to match combined channel gain
    profs['molecular'].gain_scale(mol_gain)
#    if Denoise_Mol:
#        MolDenoise.gain_scale(mol_gain)    
    
beta_aer = lp.AerosolBackscatter(profs['molecular'],profs['combined_hi'],beta_m)
    
#fig_data = lp.pcolor_profiles([beta_aer],scale=['log'],climits=[[1e-8,1e-3]]) #,ylimits=[MinAlt*1e-3,MaxAlt*1e-3])

SNR_filter = beta_aer.SNR() < 2.0
beta_aer.mask(SNR_filter)
beta_aer.mask_range('<',100)
accum_beta_a = np.nanmax(beta_aer.profile,axis=1)
plt.figure()
plt.semilogy(beta_aer.time/3600.,accum_beta_a,label=r'Max $\beta_a$')
plt.semilogy(beta_aer.time[[0,-1]]/3600.,np.array([bs_th,bs_th]),'k--',linewidth=1.5,label='Threshold')
plt.xlabel('Time [h-UTC]')
plt.ylabel('Column Aerosol Backscatter [$m^{-1}sr^{-1}$]')
plt.grid(b=True)
plt.legend()



i_rm = np.nonzero(accum_beta_a > bs_th)[0]

beta_aer.remove_time_indices(i_rm,label='Cloud Filter')
beta_m.remove_time_indices(i_rm,label='Cloud Filter')
if settings['EstimateExtinction']:
        
    alpha_aer = beta_aer.copy()
    #alpha_aer.remove_mask()
    alpha_aer.gain_scale(settings['LRassumed'],gain_var=settings['LRassumed']*0.5)
    alpha_aer.label = 'Extinction Coefficient'
    alpha_aer.descript = 'Approximated extinction coefficient based on assumed lidar ratio'
    alpha_aer.profile_type = '$m^{-1}$'
    #alpha_aer.profile[:,:150] = 0
    #alpha_aer.profile_variance[:,:150] = 0
    
    ODest = alpha_aer+8*np.pi/3.0*beta_m
    ODest.profile[:,:150] = 0
    ODest.profile_variance[:,:150] = 0
    izMax = np.argmin(np.abs(ODest.range_array-settings['ExtinctionAlt']))
    ODest.profile[:,izMax:] = 0
    ODest.profile_variance[:,izMax:]
    ODest.cumsum(axis=1)
    
    
    Tatm = np.nanmean(np.exp(-2*ODest.profile),axis=0)
    Tatm.mask= np.zeros(Tatm.shape,dtype=bool)
    Tatm_var = np.nansum(4*np.exp(-4*ODest.profile)*ODest.profile_variance,axis=0)/(ODest.time.size**2)
    Tatm_var.mask = np.zeros(Tatm_var.shape,dtype=bool)
    
    #totalExt = alpha_aer.profile+8*np.pi/3.0*beta_m.profile
    #totalExt[:,:150] = 0   # force extinction to zero below bin 110
    #
    #OD_est = sp.integrate.cumtrapz(totalExt,dx=alpha_aer.mean_dR,axis=1)
    #OD_est2 = np.nancumsum(totalExt,axis=1)*alpha_aer.mean_dR
    #Tatm0 = np.nanmean(np.exp(-2*OD_est),axis=0)
    #Tatm2 = np.nanmean(np.exp(-2*OD_est2),axis=0)
    
        
    
    plt.figure()
    plt.plot(Tatm,label='Transmission')
    plt.plot(np.sqrt(Tatm_var),label='std')
    #plt.plot(Tatm0,':',label='trapazoid')
    #plt.plot(Tatm2,'--',label='cumulative')
    plt.grid(b=True)
    plt.legend()
    #plt.plot(Tatm2,'--')
    
    Tatm = Tatm[np.newaxis,:]
else:
    Tatm = np.zeros((1,profs['molecular'].profile.shape[1]))


profs['beta_m'] = beta_m


if settings['use_bs_thr']:
    for var in profs.keys():
        if not 'Removed specified times: Cloud Filter' in profs[var].ProcessingStatus:
            profs[var].remove_time_indices(i_rm,label='Cloud Filter')
            
        profs[var].time_integrate(avg=True)
        if var != 'beta_m':
            profs[var].range_correct()
        
lp.plotprofiles(profs)

geo_raw = profs['beta_m']/profs['molecular']
geo_raw.multiply_piecewise(Tatm)

lp.plotprofiles([geo_raw],varplot=True,scale='log')


plt.show()

"""    
    

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

        save_cal_file_ez = 'diff_geo_GVHSRL'+time_dt[i0[cal_index]].strftime('%Y%m%d')+'_tmp'
        
        print('Saving to:\n  '+save_file_path+save_cal_file)        
        
        np.savetxt(save_file_path+save_cal_file,write_data,fmt='\t%d\t%.4e\t%.4e\t%.4e',header=header_str)
        
        print('Saving python vars to:\n  '+save_path_ez+save_cal_file_ez+'.npz') 
        
        np.savez(save_path_ez+save_cal_file_ez,start_time=time_dt[i0[cal_index]],stop_time=time_dt[i1[cal_index]], \
            hi_diff_geo=hi_diff_geo,lo_diff_geo=lo_diff_geo, \
            hi_diff_var=pHi.profile_variance.data.flatten(),lo_diff_var=pLo.profile_variance.data.flatten(),\
            hi_prof = pHi.profile.data.flatten(),lo_prof = pLo.profile.data.flatten(),\
            sg_win = sg_win, sg_order = sg_order, i_norm = i_norm, 
            range_array=pHi.range_array,lo_norm=lo_norm,hi_norm=hi_norm, \
            TelescopeDirection = np.nanmean(var_1d_data['TelescopeDirection'][i0[cal_index]:i1[cal_index]]))
        
"""

