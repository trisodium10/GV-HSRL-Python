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

import ExternalDataFunctions as ex
    


# input and raw_input are the same for python 3 and 2 respectively
# this makes it so input always accepts a string
try:
    input=raw_input
except NameError:
    pass



cal_file_path = os.path.abspath(__file__+'/../../calibrations/cal_files/')+'/'
cal_file = cal_file_path + 'gv_calvals.json'

reanalysis_path = os.path.abspath(__file__+'/../../../external_data/')+'/'

#save_file_path = '/Users/mhayman/Documents/Python/Lidar/'
#save_file_path = '/h/eol/mhayman/HSRL/hsrl_processing/hsrl_configuration/projDir/calfiles/'

settings = {
    'deadtime_correct':True
    }

tres = 1*60.0  # resolution in time in seconds (0.5 sec)
zres = 10.0  # resolution in altitude points (7.5 m)

#mol_gain = 1.133915#1.0728915  # gain adjustment to molecular channel

# index for where to treat the profile as background only
BGIndex = -100; # negative number provides an index from the end of the array
platform = 'ground' # 'ground' or 'airborne'.  If 'airborne' it needs an aircraft netcdf.
MaxAlt = 10e3

RemoveCals = True   # don't include instances where the I2 cell is removed
                    # scan files are not included in the file search so they
                    # are removed anyway

diff_geo_correct = True  # apply differential overlap correction

load_reanalysis = False # load T and P reanalysis from NCEP/NCAR Model

plot_2D = False   # pcolor plot the BSR and depolarization profiles

Estimate_Mol_Gain = True # use statistics on BSR to estimate the molecular gain


#sg_win = 11
#sg_order = 5

#basepath = '/scr/eldora1/HSRL_data/'  # old path - still works with link from HSRL_data to /hsrl/raw/
#basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path
basepath = '/Users/mhayman/Documents/HSRL/GVHSRL_data/'  # local computer data path
#basepath = '/scr/rain1/rsfdata/projects/socrates/hsrl/raw/' # socrates path
    
year_in = 2017
month_in = 10
day_in = 27
start_hr = 10 #17.2
stop_hr = 12# 4

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
var_2d_list = ['molecular','combined_hi','combined_lo','cross']




# grab raw data from netcdf files
[timeD,time_dt,time_sec],var_1d_data, profs = gv.load_raw_data(time_start,time_stop,var_2d_list,var_1d_list,basepath=basepath,verbose=True,as_prof=True)

# find instances in raw data where I2 cell is removed
if 'RemoveLongI2Cell' in var_1d_data.keys():
    cal_indices = np.nonzero(var_1d_data['RemoveLongI2Cell'] < 50)[0]
else:
    cal_indices = []

with open(cal_file,"r") as f:
    cal_json = json.loads(f.read())
f.close()

mol_gain,diff_geo_file = lp.get_calval(time_start,cal_json,"Molecular Gain",returnlist=['value','diff_geo'])
baseline_file = lp.get_calval(time_start,cal_json,"Baseline File")[0]
# load differential overlap correction
diff_data = np.load(cal_file_path+diff_geo_file)
baseline_data = np.load(cal_file_path+baseline_file)
#diff_data = np.load(cal_file_path+'diff_geo_GVHSRL20171025_tmp.npz')
#baseline_data = np.load(cal_file_path+'diff_geo_GVHSRL20171025_tmp.npz')
dead_time_list = lp.get_calval(time_start,cal_json,"Dead_Time",returnlist=var_2d_list)
dead_time = dict(zip(var_2d_list,dead_time_list))

lidar_location = lp.get_calval(time_start,cal_json,"Location",returnlist=['latitude','longitude'])

# set the master time to match all 2D profiles to
# (1d data will not be resampled)
master_time = np.arange(time_sec[0]-tres/2,time_sec[-1]+tres/2,tres)

time_1d,var_1d = gv.var_time_resample(master_time,time_sec,var_1d_data,average=True)
int_profs = {}  # obtain time integrated profiles
for var in profs.keys():
    if RemoveCals:
        # remove instances where the I2 cell is removed
        profs[var].remove_time_indices(cal_indices)
    if settings['deadtime_correct']:
        if hasattr(profs[var],'NumProfsList'):
            profs[var].nonlinear_correct(dead_time[var],laser_shot_count=2000*profs[var].NumProfsList[:,np.newaxis],std_deadtime=5e-9)
        else:
            # number of laser shots is based on an assumption that there is one 0.5 second profile per time bin
            profs[var].nonlinear_correct(dead_time[var],laser_shot_count=2000,std_deadtime=5e-9)
    profs[var].time_resample(tedges=master_time,update=True,remainder=False)
    int_profs[var] = profs[var].copy()
    int_profs[var].time_integrate()
    
    profs[var].bg_subtract(BGIndex)
    if var == 'combined_hi' and diff_geo_correct:
        profs[var].diff_geo_overlap_correct(diff_data['hi_diff_geo'])
    elif var == 'combined_lo' and diff_geo_correct:
        profs[var].diff_geo_overlap_correct(diff_data['lo_diff_geo'])
#        profs[var].gain_scale(1.0/diff_data['lo_norm'])
    profs[var].slice_range(range_lim=[0,MaxAlt])
    int_profs[var].bg_subtract(BGIndex)
    int_profs[var].slice_range(range_lim=[0,MaxAlt])


profs['combined'],_ = gv.merge_hi_lo(profs['combined_hi'],profs['combined_lo'],plot_res=True)

if load_reanalysis:
    pres,temp = ex.load_fixed_point_NCEP_TandP(profs['molecular'],lidar_location,reanalysis_path)


lp.plotprofiles(profs)

profs['molecular'].gain_scale(mol_gain)
#profs['combined_hi'].diff_geo_overlap_correct(diff_data['hi_diff_geo'][:profs['combined_hi'].range_array.size])
#profs['combined_lo'].diff_geo_overlap_correct(diff_data['lo_diff_geo'][:profs['combined_lo'].range_array.size])
#profs['combined_lo'].gain_scale(1.0/diff_data['lo_norm'])

BSR = profs['combined_hi']/profs['molecular']
BSR.descript = 'Ratio of combined to molecular backscatter'
BSR.label = 'Backscatter Ratio'
BSR.profile_type = 'unitless'

#BSR_mask = (BSR.profile-1)/np.sqrt(BSR.profile_variance) < 5.0
BSR_mask = BSR.profile < 2.0

#BSR2 = profs['combined_hi'].copy()
#BSR2.divide_prof(profs['molecular'])
#BSR2.descript = 'Ratio of combined to molecular backscatter'
#BSR2.label = 'Backscatter Ratio'
#BSR2.profile_type = 'unitless'

dVol = profs['cross']/(profs['combined_hi']+profs['cross'])
#dVol = profs['combined_hi'].copy()
dVol.descript = 'Propensity of Volume to depolarize (d)'
dVol.label = 'Volume Depolarization'
dVol.profile_type = 'unitless'

d_mol = 2*0.000365/(1+0.000365) # molecular depolarization

#Particle Depolarization = dVol/(1.0-1.0/BSR) - d_mol/(BSR-1)
dPart = (BSR*dVol-d_mol)/(BSR-1)
dPart.descript = 'Propensity of Particles to depolarize (d)'
dPart.label = 'Particle Depolarization'
dPart.profile_type = 'unitless'


if Estimate_Mol_Gain:
    # This segment estimates what the molecular gain should be 
    # based on a histogram minimum in BSR over the loaded data
    
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


# add a diagnostic for counts/backscatter coeff

# add a diagnostic for diff overlap between lo and hi channels as a function
# of count rate or backscatter coeff

dPartMask = dPart.SNR() < 3.0


if plot_2D:
    lp.pcolor_profiles([BSR,dPart],scale=['log','linear'],climits=[[1,5e2],[0,0.7]])
    #lp.plotprofiles(profs)
    #dPart.mask(dPartMask)
    #lp.pcolor_profiles([BSR,dVol],scale=['log','linear'],climits=[[1,5e2],[0,1.0]])
    #lp.pcolor_profiles([dVol],scale=['linear'],climits=[[0,1]])

#import scipy.optimize
#errfun = lambda x: np.nansum((profs['combined_hi'].profile.flatten()-profs['combined_lo'].profile.flatten()*x)**2/(profs['combined_hi'].profile_variance.flatten()+profs['combined_lo'].profile_variance.flatten()*x**2))
#x0 =  1.0
#sol1D = scipy.optimize.minimize_scalar(errfun) #,disp=0,x0,maxfun=2000,eta=1e-5 , ,opt_iterations,opt_exit_mode
#comb_lo_gain = sol1D['x']
#profs['combined_lo'].gain_scale(comb_lo_gain)
#
#
#plt.figure(); 
#plt.scatter(np.log10(profs['combined_lo'].profile.flatten()),np.log10(profs['combined_hi'].profile.flatten())); 
#plt.plot([np.nanmin(np.log10(profs['combined_lo'].profile.flatten())),np.nanmax(np.log10(profs['combined_lo'].profile.flatten()))],[np.nanmin(np.log10(profs['combined_lo'].profile.flatten())),np.nanmax(np.log10(profs['combined_lo'].profile.flatten()))],'k--')
#
#combined = profs['combined_hi'].copy()
#combined_to_lo = np.nonzero(np.logical_or( 
#    profs['combined_hi'].profile_variance > profs['combined_lo'].profile_variance,
#    np.isnan(profs['combined_hi'].profile)))
#combined.profile[combined_to_lo]= profs['combined_lo'].profile[combined_to_lo]
#combined.profile_variance[combined_to_lo]= profs['combined_lo'].profile_variance[combined_to_lo]
#combined.descript = 'Merged hi/lo gain combined channel'
#combined.label = 'Merged Combined Channel'
#
#lp.plotprofiles([profs['combined_lo'],profs['combined_hi'],combined],time=18.4*3600.0,varplot=True)