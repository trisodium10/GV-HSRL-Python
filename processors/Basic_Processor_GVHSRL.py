# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 15:15:04 2017

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
import LidarProfileFunctions as lp
import datetime
#import glob

import json

import GVHSRLlib as gv
    


# input and raw_input are the same for python 3 and 2 respectively
# this makes it so input always accepts a string
try:
    input=raw_input
except NameError:
    pass



cal_file_path = os.path.abspath(__file__+'/../../calibrations/cal_files/')+'/'
cal_file = cal_file_path + 'gv_calvals.json'
#save_file_path = '/Users/mhayman/Documents/Python/Lidar/'
#save_file_path = '/h/eol/mhayman/HSRL/hsrl_processing/hsrl_configuration/projDir/calfiles/'

tres = 1*60.0  # resolution in time in seconds (0.5 sec)
zres = 10.0  # resolution in altitude points (7.5 m)

#mol_gain = 1.133915#1.0728915  # gain adjustment to molecular channel

# index for where to treat the profile as background only
BGIndex = -100; # negative number provides an index from the end of the array
platform = 'ground' # 'ground' or 'airborne'.  If 'airborne' it needs an aircraft netcdf.
MaxAlt = 14e3

RemoveCals = True   # don't include instances where the I2 cell is removed
                    # scan files are not included in the file search so they
                    # are removed anyway

diff_geo_correct = True  # apply differential overlap correction

sg_win = 11
sg_order = 5

#basepath = '/scr/eldora1/HSRL_data/'  # old path - still works with link from HSRL_data to /hsrl/raw/
basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path
#basepath = '/Users/mhayman/Documents/HSRL/GVHSRL_data/'  # local computer data path
        
year_in = 2017
month_in = 10
day_in = 27
start_hr = 17.2
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

time_start = datetime.datetime(year_in,month_in,day_in)+datetime.timedelta(hours=start_hr)
time_stop = time_start + datetime.timedelta(hours=stop_hr)


# list of 1D variables to load
var_1d_list = ['total_energy','RemoveLongI2Cell'\
    ,'TelescopeDirection','TelescopeLocked','polarization','DATA_shot_count']  # 'DATA_shot_count'

# list of 2D variables (profiles) to load
var_2d_list = ['molecular','combined_hi','combined_lo']  # 'cross',






# grab raw data from netcdf files
[timeD,time_dt,time_sec],var_1d_data, profs = gv.load_raw_data(time_start,time_stop,var_2d_list,var_1d_list,basepath=basepath,verbose=True,as_prof=True)

# find instances in raw data where I2 cell is removed
if 'RemoveLongI2Cell' in var_1d_data.keys():
    cal_indices = np.nonzero(var_1d_data['RemoveLongI2Cell'] < 50)[0]
else:
    cal_indices = []

with open(cal_file,"r") as f:
    cal_json = json.loads(f.read())

mol_gain,diff_geo_file = lp.get_calval(time_start,cal_json,"Molecular Gain",returnlist=['value','diff_geo'])
baseline_file = lp.get_calval(time_start,cal_json,"Baseline File")[0]
# load differential overlap correction
diff_data = np.load(cal_file_path+diff_geo_file)
baseline_data = np.load(cal_file_path+baseline_file)
#diff_data = np.load(cal_file_path+'diff_geo_GVHSRL20171025_tmp.npz')
#baseline_data = np.load(cal_file_path+'diff_geo_GVHSRL20171025_tmp.npz')


# set the master time to match all 2D profiles to
# (1d data will not be resampled)
master_time = np.arange(time_sec[0]-tres/2,time_sec[-1]+tres/2,tres)

time_1d,var_1d = gv.var_time_resample(master_time,time_sec,var_1d_data,average=True)
int_profs = {}  # obtain time integrated profiles
for var in profs.keys():
    if RemoveCals:
        # remove instances where the I2 cell is removed
        profs[var].remove_time_indices(cal_indices)
    profs[var].time_resample(tedges=master_time,update=True,remainder=False)
    int_profs[var] = profs[var].copy()
    int_profs[var].time_integrate()
    
    profs[var].bg_subtract(BGIndex)
    if var == 'combined_hi' and diff_geo_correct:
        profs[var].diff_geo_overlap_correct(diff_data['hi_diff_geo'])
    elif var == 'combined_lo' and diff_geo_correct:
        profs[var].diff_geo_overlap_correct(diff_data['hi_diff_geo'])
        profs[var].gain_scale(1.0/diff_data['lo_norm'])
    profs[var].slice_range(range_lim=[0,MaxAlt])
    int_profs[var].bg_subtract(BGIndex)
    int_profs[var].slice_range(range_lim=[0,MaxAlt])







lp.plotprofiles(profs)

profs['molecular'].gain_scale(mol_gain)
#profs['combined_hi'].diff_geo_overlap_correct(diff_data['hi_diff_geo'][:profs['combined_hi'].range_array.size])
#profs['combined_lo'].diff_geo_overlap_correct(diff_data['lo_diff_geo'][:profs['combined_lo'].range_array.size])
#profs['combined_lo'].gain_scale(1.0/diff_data['lo_norm'])

BSR = profs['combined_hi'].copy()
BSR.descript = 'Ratio of combined to molecular backscatter'
BSR.label = 'Backscatter Ratio'
BSR.profile_type = 'unitless'
BSR.divide_prof(profs['molecular'])

bbsr = np.linspace(0,4,400)
bsnr = np.linspace(10,150,100)
hbsr = np.histogram2d(BSR.profile.flatten(),BSR.SNR().flatten(),bins=[bbsr,bsnr])
#plt.figure()
#plt.pcolor(bbsr,bsnr,hbsr[0].T)

i_hist_median = np.argmax(hbsr[0],axis=0)
iset = np.arange(hbsr[0].shape[1])
dh1 = hbsr[0][i_hist_median,iset]-hbsr[0][i_hist_median-1,iset]
dh2 = hbsr[0][i_hist_median+1,iset]-hbsr[0][i_hist_median,iset]
dbsr = np.mean(np.diff(bbsr))
bsr1 = bbsr[i_hist_median]
bsr2 = bbsr[i_hist_median+1]

m_0 = (dh2-dh1)/dbsr
b_0 = dh1-m_0*bsr1

Nsm = 6  # number of bins to smooth over
hist_median = (-b_0)/m_0
hist_med_sm = np.convolve(hist_median,np.ones(Nsm)*1.0/Nsm,mode='same')

#hist_median = bbsr[i_hist_median]

plt.figure()
plt.pcolor(bbsr,bsnr,hbsr[0].T)
plt.plot(hist_median,bsnr[1:],'r--')
plt.plot(hist_med_sm,bsnr[1:],'g--')
plt.xlabel('BSR')
plt.ylabel('BSR SNR')

i_snr_lim = np.nonzero(np.logical_and(bsnr < 80,bsnr > 20))[0]

mol_gain_adj = np.nanmin(hist_med_sm[i_snr_lim])

print('Current Molecular Gain: %f'%mol_gain)
print('Suggested Molecular Gain: %f'%(mol_gain*mol_gain_adj))

#lp.plotprofiles(profs)
lp.pcolor_profiles([BSR],scale=['log'],climits=[[1,5e2]])
