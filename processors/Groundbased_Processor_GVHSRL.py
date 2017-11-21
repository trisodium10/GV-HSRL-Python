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

import MLELidarProfileFunctions as mle

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

load_reanalysis = True # load T and P reanalysis from NCEP/NCAR Model

plot_2D = True   # pcolor plot the BSR and depolarization profiles

Estimate_Mol_Gain = True # use statistics on BSR to estimate the molecular gain

Denoise_Mol = False  # run PTV denoising on molecular channel

hsrl_rb_adjust = True # apply rayleigh brillouin correction to molecular profile


#sg_win = 11
#sg_order = 5

#basepath = '/scr/eldora1/HSRL_data/'  # old path - still works with link from HSRL_data to /hsrl/raw/
basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path
#basepath = '/Users/mhayman/Documents/HSRL/GVHSRL_data/'  # local computer data path
        
year_in = 2017
month_in = 10
day_in = 27
start_hr = 16 #17.2
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
if hsrl_rb_adjust:
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
    
    if var == 'molecular':
        MolRaw = profs['molecular'].copy()
    
    profs[var].bg_subtract(BGIndex)
    if var == 'combined_hi' and diff_geo_correct:
        profs[var].diff_geo_overlap_correct(diff_data['hi_diff_geo'])
    elif var == 'combined_lo' and diff_geo_correct:
        profs[var].diff_geo_overlap_correct(diff_data['hi_diff_geo'])
        profs[var].gain_scale(1.0/diff_data['lo_norm'])
    
    profs[var].slice_range(range_lim=[0,MaxAlt])
    int_profs[var].bg_subtract(BGIndex)
    int_profs[var].slice_range(range_lim=[0,MaxAlt])




if load_reanalysis:
    pres,temp = ex.load_fixed_point_NCEP_TandP(profs['molecular'],lidar_location,reanalysis_path)
    beta_m = lp.get_beta_m(temp,pres,profs['molecular'].wavelength)

if Denoise_Mol:
#    tune_data = mle.DenoiseBG(MolRaw,-10,verbose=False,plot_sol=True,tv_lim =[0.4, 1.4],N_tv_pts=78)
    MolDenoise,tune_list = mle.DenoiseMolecular(MolRaw,beta_m_sonde=beta_m, \
                            MaxAlt=MaxAlt,accel = False,tv_lim =[1.5, 2.8],N_tv_pts=59, \
                            geo_data=dict(geo_prof=np.array([2e14])),bg_index=-10)
    


if hsrl_rb_adjust:
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
    
    if Denoise_Mol:
        MolDenoise.multiply_piecewise(1.0/eta_i2)
        MolDenoise.gain_scale(mol_gain)
else:
    # Rescale molecular channel to match combined channel gain
    profs['molecular'].gain_scale(mol_gain)
    if Denoise_Mol:
        MolDenoise.gain_scale(mol_gain)


lp.plotprofiles(profs)

#profs['molecular'].gain_scale(mol_gain)
#
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

beta_aer = lp.AerosolBackscatter(profs['molecular'],profs['combined_hi'],beta_m)
if Denoise_Mol:
    beta_aer_denoise = lp.AerosolBackscatter(MolDenoise,profs['combined_hi'],beta_m)

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

#dPartMask = dPart.SNR() < 3.0


if plot_2D:
    lp.pcolor_profiles([beta_aer,dPart],scale=['log','linear'],climits=[[1e-8,1e-4],[0,0.7]])
    if Denoise_Mol:
        lp.pcolor_profiles([beta_aer_denoise],scale=['log'],climits=[[1e-8,1e-4]])
#    lp.pcolor_profiles([BSR,dPart],scale=['log','linear'],climits=[[1,5e2],[0,0.7]])
    #lp.plotprofiles(profs)
    #dPart.mask(dPartMask)
    #lp.pcolor_profiles([BSR,dVol],scale=['log','linear'],climits=[[1,5e2],[0,1.0]])
    #lp.pcolor_profiles([dVol],scale=['linear'],climits=[[0,1]])

plt.show()