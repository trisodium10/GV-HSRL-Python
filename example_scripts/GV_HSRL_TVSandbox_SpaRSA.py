# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 08:20:36 2018

@author: mhayman
"""

import numpy as np
#import matplotlib.pyplot as plt
import scipy.optimize
#from mpl_toolkits.mplot3d import Axes3D

import LidarProfileFunctions as lp
import LidarPlotFunctions as lplt
import GVHSRLlib as gv
import MLELidarProfileFunctions as mle

import datetime

import netCDF4 as nc4

import json

import matplotlib

import os

#matplotlib.use('Agg')
    
import matplotlib.pyplot as plt    
from mpl_toolkits.mplot3d import Axes3D 
   
# path to files being analyzed
#ncpath = '/scr/sci/mhayman/SOCRATES/range_centered/test_data/SOCRATESrf03/'
##ncpath = '/Users/mhayman/Documents/HSRL/GVHSRL/test_data/SOCRATESrf03/'
## list of files to analyze
   #2DC cloud profiles
#nclist0 = ['SOCRATESrf03_GVHSRL_20180123T0210_20180123T0220.nc']
#t_lim = [26*3600+12*60,26*3600+14*60] #
#r_lim = [0,2e3]
#step_eps= 1e-5#2e-5
#max_iter = 1000
#lam_set = {'xB':69.9,'xS':6.26e-3,'xP':4.32}
#lam_range = {'xB':[1.5,2.5],'xS':[-2.5,1.5],'xP':[0,1.0]}
#Num_reg_iter = 1


# oriented ice crystal data
ncpath = '/scr/sci/mhayman/SOCRATES/range_centered/test_data/SOCRATESrf07/'
nclist0 = ['SOCRATESrf07_GVHSRL_20180131T0610_20180131T0630.nc']
t_lim = [6.42*3600,6.44*3600] #
r_lim = [1e3,5.5e3]
step_eps= 2e-9#2e-5
max_iter = 10000
lam_set = {'xB':85.0,'xS':0.2,'xP':55.11} #{'xB':85.0,'xS':0.2,'xP':1.11}
lam_range = {'xB':[1.0,2.5],'xS':[-2.5,0],'xP':[0,1.5]}
Num_reg_iter = 1

# if desired, loaded files will be constrainted by time
time_start = datetime.datetime(year=1900,month=1,day=1)
time_stop = datetime.datetime(year=2100,month=1,day=1)

load_mask = False
save_results = True

#t_lim = [26*3600+12*60,26*3600+14*60] #
#r_lim = [0,2e3]
#step_eps= 1e-6#2e-5
#max_iter = 1000



#savefigpath = '/h/eol/mhayman/HSRL/Plots/PTV/20180118_26UTC_20180519/'
#
#savencfile='20180118_26UTC_20180519.nc'

# set the variable bounds in their native space.
# this is converted to the state variable later
opt_bounds = {'xB':[1e-9,1e-2],
              'xS':[1.0,100.0],
              'xP':[0.0,1.0]}

#prof_list = ['Aerosol_Backscatter_Coefficient','Aerosol_Extinction_Coefficient','Merged_Combined_Channel',
#             'Particle_Depolarization','Volume_Depolarization','Denoised_Aerosol_Backscatter_Coefficient'] #,'Denoised_Aerosol_Backscatter_Coefficient'] #'Aerosol_Extinction_Coefficient'

prof_list = ['Aerosol_Backscatter_Coefficient','Particle_Depolarization','Aerosol_Extinction_Coefficient',
             'Raw_Molecular_Backscatter_Channel','Raw_High_Gain_Total_Backscatter_Channel',
             'Raw_Low_Gain_Total_Backscatter_Channel','Raw_Cross_Polarization_Channel',
             'Temperature','Pressure','Molecular_Backscatter_Coefficient']

air_vars = ['GGALT','ROLL','PITCH','GGLAT','GGLON']  # aircraft variables to plot on the figures

lidar_vars = ['lidar_pointing','polarization','TelescopeDirection']


plot_settings = {
    'Aerosol_Backscatter_Coefficient':{
        'climits':[1e-8,1e-3],
        'scale':'log',
	'colormap':'jet'
        },
    'Aerosol_Extinction_Coefficient':{
        'climits':[1e-5,1e-2],
        'scale':'log',
	'colormap':'viridis'
        },
    'Particle_Depolarization':{
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
        'climits':[1e-8,1e-3],
        'scale':'log',
	'colormap':'viridis'
        },
    'Merged_Combined_Channel':{
        'climits':[1e-1,1e4],
        'scale':'log',
	'colormap':'jet'
        },
    'Lidar_Ratio':{
        'climits':[15,40],
        'scale':'linear',
	'colormap':'viridis'
        },
    }


#%%  Load Data Files

"""
Load selected data files
"""

# sort the files by name and add the path to the files
nclist = [ncpath+x for x in sorted(nclist0)]


findex = []
aircraft_data = {}
lidar_data = {}

# Code to constrain the files by time.
for ai in range(len(nclist)):
    filestart = datetime.datetime.strptime(nclist[ai][-30:-17],'%Y%m%dT%H%M')
    filestop = datetime.datetime.strptime(nclist[ai][-16:-3],'%Y%m%dT%H%M')

    if time_start <= filestart and time_stop >= filestart:
        findex.extend([ai])
    elif time_start <= filestop and time_stop >= filestart:
        findex.extend([ai])
    

#save_plots_base = process_vars['flt']+'_GVHSRL_'+time_start.strftime('%Y%m%dT%H%M')+'_'+time_stop.strftime('%Y%m%dT%H%M')    
#tlims = [(time_start-day_start).total_seconds()/3600.0,(time_stop-day_start).total_seconds()/3600.0]
#print(time_start.strftime('%H:%M') + ' - ' + time_stop.strftime('%H:%M')+':')


prof_start = True
profs = {}
aircraft_data = {}
lidar_data = {}
for fi in findex:
    ncfilename = nclist[fi]
    f = nc4.Dataset(ncfilename,'r')     
    for a_var in air_vars:
#            print(a_var)
        if a_var in aircraft_data.keys():
#                print('... concatenating')
            aircraft_data[a_var] = np.concatenate((aircraft_data[a_var],lp.ncvar(f,a_var)))
        else:
#                print('... first find')
            aircraft_data[a_var] = lp.ncvar(f,a_var)
    
    # load lidar variables
    for l_var in lidar_vars:
        new_data = lp.ncvar(f,l_var)
        if l_var in lidar_data.keys():                
            if new_data.ndim == 2:
                lidar_data[l_var] = np.hstack((lidar_data[l_var],new_data))
            else:
                lidar_data[l_var] = np.concatenate((lidar_data[l_var],new_data))
        else:
            lidar_data[l_var] = new_data
            

    t_data0 = lp.ncvar(f,'time').astype(np.float)
    f.close()
    

    
    
    if prof_start:
        t_data = t_data0.copy()
        prof_start = False
    else:
        t_data = np.concatenate((t_data,t_data0))
    
    for var in prof_list:
        prof0 = lp.load_nc_Profile(ncfilename,var,mask=load_mask)
	    # check if the load was successful
        if hasattr(prof0,'time'):
            if var in profs.keys():
                profs[var].cat_time(prof0,front=False)
            else:
                profs[var] = prof0.copy()
                
#%% configure save path and file names
var = prof_list[0]
if save_results:
    save_path = profs[var].StartDate.strftime('/h/eol/mhayman/HSRL/Plots/PTV/Denoise_%Y%m%d/')
    if not os.path.exists(save_path):
                os.makedirs(save_path)
    save_file = profs[var].StartDate.strftime('%Y%m%d')+'_%d_%d_UTC_'%(t_lim[0]/3600,t_lim[1]/3600)+datetime.datetime.now().strftime('%Y%m%d_%H%M')
    savefigpath = save_path + save_file
    savencfile = save_file +'.nc'
#%%  Plot Sample Profile

#var = 'Aerosol_Backscatter_Coefficient'
#lplt.scatter_z(profs[var],lidar_pointing=lidar_data['lidar_pointing'],
#               lidar_alt=aircraft_data['GGALT'],climits=plot_settings[var]['climits'],
#               cmap=plot_settings[var]['colormap'],scale=plot_settings[var]['scale'],s=2,t_axis_scale=5.0)

#%%  Trim and prep data for processing

"""
Prep data for processing
"""

#def MLE_Cals_2D(MolRaw,CombRaw,beta_aer,surf_temp,surf_pres,geo_data,minSNR=0.5,\
#t_lim=np.array([np.nan,np.nan]),verify=False,verbose=False,print_sol=True,\
#plotfigs=True,lam_array=np.array([np.nan]),Nmax=2000):
"""
Runs a maximum likelihood estimator to correct for
    Channel gain mismatch, 
    Aerosol to Molecular Channel Crosstalk
    Detector Dead Time
and obtain estimates of
    Backscatter coefficient
    Extinction coefficient
    Lidar Ratio
    Each Channel's gain
    Aerosol to Molecular Channel Crosstalk
    Detector Dead Time
    
Inputs:
    MolRaw - Raw Molecular Profile
    CombRaw - Raw Combined Profile
    beta_aer - derived estimate of aerosol backscatter coefficent
    surf_temp - temperature data from the surface station
    surf_pres - pressure data from the surface station
    geo_data - data file from loading the geometric overlap function
    minSNR - minimum SNR required to be included as containing possible aerosols
    t_lim - currently disabled.  when added, will allow us to select only a segment to operate on.
    verify - if True it will use Poisson Thinning to verify TV solution
    verbose - if True it will output the fit result of each TV iteration
    print_sol - print the solved calibration parameters
    plotfigs - plot the TV verification error
    lam_array - array containing TV values to be evaluated
    Nmax - maximum optimizer iterations
    
Substitutions:
    beta_aer for aer_beta_dlb
    beta_aer_E for aer_beta_E
"""

# Trim the data in time and range

var = 'Raw_Molecular_Backscatter_Channel'
t_data = profs[var].time.copy()  # store time data for 1d variables

raw_range = profs[var].range_array.copy()  # used to trim the geo files later

#t_lim = [26*3600+12*60,26*3600+14*60] #
#r_lim = [0,2e3]

# force range limits to as or more limiting than most restrictive profile
for var in profs.keys():
    if r_lim[0] < profs[var].range_array[0]:
        r_lim[0] = profs[var].range_array[0]
    if r_lim[1] > profs[var].range_array[-1]:
        r_lim[1] = profs[var].range_array[-1]

# additional processing for the raw profiles:
# 1.   Estimate background, but don't subtract it
# 2.   Align the range grid with process profiles
for var in profs.keys():
    profs[var].slice_time(t_lim)
    if 'Raw' in var:
        # run a custom background estimation
        # store the result in the profile background (bg) but don't
        # actually subtract it
        profs[var].bg = np.nanmean(profs[var].profile[:,-200:],axis=1)
        # now trim the range dimension to match the processed profiles
#        profs[var].slice_range(range_lims=[profs['Aerosol_Backscatter_Coefficeint'].range_array[0],profs['Aerosol_Backscatter_Coefficeint'].range_array[1]])
    profs[var].slice_range(range_lim=r_lim)


# trim time on 1D variables
irm1d = np.nonzero((t_data < t_lim[0]) +(t_data > t_lim[1]))
for var1d in aircraft_data.keys():
    aircraft_data[var1d] = np.delete(aircraft_data[var1d],irm1d)
for var1d in lidar_data.keys():
    if lidar_data[var1d].ndim > 1:
        lidar_data[var1d] = np.delete(lidar_data[var1d],irm1d,axis=1)
    else:
        lidar_data[var1d] = np.delete(lidar_data[var1d],irm1d)
    

# assign an x,y,z and t value for each lidar pixel

var = 'Raw_Molecular_Backscatter_Channel'
# time of each pixel
l_time = (np.ones((1,profs[var].range_array.size))*(profs[var].time[:,np.newaxis]))

# get aircraft (x,y) from latitude and longitude data
a_pos_x,a_pos_y = lp.dist_from_latlon(aircraft_data['GGLAT'],aircraft_data['GGLON'])
a_pos_z = aircraft_data['GGALT'].copy()

# get (x,y,z) position of every lidar pixel
l_pos_x = (a_pos_x[:,np.newaxis]+lidar_data['lidar_pointing'][0,:][:,np.newaxis]*profs[var].range_array)
l_pos_y = (a_pos_y[:,np.newaxis]+lidar_data['lidar_pointing'][1,:][:,np.newaxis]*profs[var].range_array)
l_pos_z = (a_pos_z[:,np.newaxis]-lidar_data['lidar_pointing'][2,:][:,np.newaxis]*profs[var].range_array)

# copy the raw profiles to avoid modifying them
raw_profs = {}
for var in profs.keys():
    if 'Raw' in var:
        raw_profs[var] = profs[var].copy()

#%%  Load cals and trim them
"""
Load Calibration Parameters
"""

cal_file_path = '/h/eol/mhayman/PythonScripts/HSRL_Processing/GV-HSRL-Python/calibrations/cal_files/'
#cal_file_path = '/Users/mhayman/Documents/Python/Lidar/GV-HSRL-Python/calibrations/cal_files/'
cal_file = cal_file_path+'gv_calvals.json'

with open(cal_file,"r") as f:
    cal_json = json.loads(f.read())
f.close()

# This isn't that sophisticated, but it works.  Don't try to update the time_start variable
time_start = profs['Raw_Molecular_Backscatter_Channel'].StartDate

# Load the appropriate calibrations
mol_gain_up,diff_geo_file,mol_gain_down,diff_geo_file_down = lp.get_calval(time_start,cal_json,'Molecular Gain',cond=[['RB_Corrected','=','True']],returnlist=['value','diff_geo','down_gain','down_diff'])  
baseline_file = lp.get_calval(time_start,cal_json,"Baseline File")[0]
diff_pol_file,pol_cal_file = lp.get_calval(time_start,cal_json,"Polarization",returnlist=['diff_geo','cal_file'])
i2_file = lp.get_calval(time_start,cal_json,"I2 Scan")
dead_time_list = lp.get_calval(time_start,cal_json,"Dead_Time",returnlist=['combined_hi','cross','combined_lo','molecular'])
dead_time = dict(zip(['High_Gain_Total_Backscatter_Channel','Cross_Polarization_Channel','Low_Gain_Total_Backscatter_Channel','Molecular_Backscatter_Channel'],dead_time_list))


geo_file_up,geo_file_down = lp.get_calval(time_start,cal_json,"Geo File",returnlist=['value','down_file'])
geo_up = np.load(cal_file_path+geo_file_up)
if len(geo_file_down) > 0:
    geo_down = np.load(cal_file_path+geo_file_down)
else:
    geo_data = geo_up
                       
# load differential overlap correction
diff_data_up = np.load(cal_file_path+diff_geo_file)
if len(diff_geo_file_down) > 0:
    diff_data_down = np.load(cal_file_path+diff_geo_file_down)
else:
    diff_data = diff_data_up
    
if len(diff_pol_file):
    diff_pol_data = np.load(cal_file_path+diff_pol_file)

diff_pol_geo = np.interp(profs['Raw_Molecular_Backscatter_Channel'].range_array,raw_range,diff_pol_data['cross_diff_geo'])


# load i2 scan from file
if len(i2_file):
    i2_data = np.load(cal_file_path+i2_file[0])
    
with open(cal_file_path+pol_cal_file,"r") as f:
    pol_cals = json.loads(f.read())
f.close()

#pol_cals['calibration']['qwp_rot_offset'] = 0  # set the offset to zero.  This is for when its rotating

# Get the transmitted Stokes vector and each channel's detected Diattenuation vectors
[Stx,DxV,DtV,DmV,DlV]=gv.HSRLPolConfig(pol_cals['calibration'],np.array([np.mean(lidar_data['polarization'])]))

#DxV = np.array([[1,0,0,-1]])
# Elements needed to compute polarization dependence of each channel
# backscatter = beta*(pol[var][0]+p*pol[var][1]) where p = 1-d
pol = {}
pol['Raw_Low_Gain_Total_Backscatter_Channel'] = gv.HSRL_FormMeasMatrix(Stx,DlV,2)
pol['Raw_High_Gain_Total_Backscatter_Channel'] = gv.HSRL_FormMeasMatrix(Stx,DtV,2)
pol['Raw_Molecular_Backscatter_Channel'] = gv.HSRL_FormMeasMatrix(Stx,DmV,2)
pol['Raw_Cross_Polarization_Channel'] = gv.HSRL_FormMeasMatrix(Stx,DxV,2)

dmol = 2*0.00365/(1+0.00365)
pmol = 1-dmol


# setup variable diffierential overlap (up vs down pointing) if supplied
if len(diff_geo_file_down) > 0:
    diff_data = {}
    key_list = ['hi_diff_geo','lo_diff_geo']
    for var in key_list:
        diff_data[var] = np.ones((lidar_data['TelescopeDirection'].size,profs['Raw_Molecular_Backscatter_Channel'].range_array.size))
        diff_data[var][np.nonzero(lidar_data['TelescopeDirection']==1.0)[0],:] = np.interp(profs['Raw_Molecular_Backscatter_Channel'].range_array,raw_range,diff_data_up[var])
        diff_data[var][np.nonzero(lidar_data['TelescopeDirection']==0.0)[0],:] = np.interp(profs['Raw_Molecular_Backscatter_Channel'].range_array,raw_range,diff_data_down[var])
   

# setup variable geo overlap (up vs down pointing) if supplied
if len(geo_file_down) > 0:
    geo_data = {}
    key_list = ['geo_mol','geo_mol_var','Nprof']
    for var in key_list:
        if var in geo_up.keys():
            geo_data[var] = np.ones((lidar_data['TelescopeDirection'].size,profs['Raw_Molecular_Backscatter_Channel'].range_array.size))
            geo_data[var][np.nonzero(lidar_data['TelescopeDirection']==1.0)[0],:] = np.interp(profs['Raw_Molecular_Backscatter_Channel'].range_array,raw_range,geo_up[var])
            if var in geo_down.keys():
                geo_data[var][np.nonzero(lidar_data['TelescopeDirection']==0.0)[0],:] = np.interp(profs['Raw_Molecular_Backscatter_Channel'].range_array,raw_range,geo_down[var])
            else:
                geo_data[var][np.nonzero(lidar_data['TelescopeDirection']==0.0)[0],:] = np.interp(profs['Raw_Molecular_Backscatter_Channel'].range_array,raw_range,geo_up[var])
        else:
            geo_data[var] = np.ones((lidar_data['TelescopeDirection'].size,1))

print('Obtaining Rayleigh-Brillouin Correction')
dnu = 20e6  # resolution
nu_max = 10e9 # max frequency relative to line center
nu = np.arange(-nu_max,nu_max,dnu)
Ti2 = np.interp(nu,i2_data['freq']*1e9,i2_data['mol_scan'])  # molecular transmission
Tam = np.interp(0,i2_data['freq']*1e9,i2_data['mol_scan'])  # aersol transmission into molecular channel

Tc2 = np.interp(nu,i2_data['freq']*1e9,i2_data['combined_scan'])  # combined transmission
Tac = np.interp(0,i2_data['freq']*1e9,i2_data['combined_scan'])  # aersol transmission into combined channel

[eta_i2,eta_c] = lp.RB_Efficiency([Ti2,Tc2],profs['Temperature'].profile.flatten(),profs['Pressure'].profile.flatten()*9.86923e-6,profs['Raw_Molecular_Backscatter_Channel'].wavelength,nu=nu,norm=True,max_size=10000)

eta_i2 = eta_i2.reshape(profs['Temperature'].profile.shape)
eta_c = eta_c.reshape(profs['Pressure'].profile.shape)

# Transmission through the atmosphere due to molecular scattering
Tatm_mol = np.exp(-2*np.cumsum(8*np.pi/3*profs['Molecular_Backscatter_Coefficient'].profile,axis=1)*profs['Molecular_Backscatter_Coefficient'].mean_dR)

#%% Build dictionary of constant range multipliers
ConstTerms0 = {}
for var in raw_profs.keys():
    ConstTerms0[var] = {'mult':np.zeros(raw_profs[var].profile.shape),'mol':np.zeros(raw_profs[var].profile.shape),'bg':np.zeros((raw_profs[var].time.size,1)),'pol':np.zeros(2),'Ta':1.0}
    ConstTerms0[var]['mult'] = Tatm_mol/geo_data['geo_mol']/raw_profs[var].range_array[np.newaxis,:]**2
    ConstTerms0[var]['bg'] = raw_profs[var].bg[:,np.newaxis]
    ConstTerms0[var]['pol'] = pol[var]
    if 'High' in var:
        ConstTerms0[var]['mult'] = Tac*ConstTerms0[var]['mult']/diff_data['hi_diff_geo']
        ConstTerms0[var]['mol'] = eta_c/Tac*profs['Molecular_Backscatter_Coefficient'].profile*(pol[var][0]+pol[var][1]*pmol)
        ConstTerms0[var]['Ta'] = Tac
    elif 'Low' in var:
        ConstTerms0[var]['mult'] = Tac*ConstTerms0[var]['mult']/diff_data['lo_diff_geo']
        ConstTerms0[var]['mol'] = eta_c/Tac*profs['Molecular_Backscatter_Coefficient'].profile*(pol[var][0]+pol[var][1]*pmol)
        ConstTerms0[var]['Ta'] = Tac
    elif 'Cross' in var:
        ConstTerms0[var]['mult'] = Tac*ConstTerms0[var]['mult']/diff_data['hi_diff_geo']/diff_pol_geo[np.newaxis,:]
        ConstTerms0[var]['mol'] = eta_c/Tac*profs['Molecular_Backscatter_Coefficient'].profile*(pol[var][0]+pol[var][1]*pmol)
        ConstTerms0[var]['Ta'] = Tac
    elif 'Molecular' in var:
        ConstTerms0[var]['mult'] = Tam*ConstTerms0[var]['mult']
        ConstTerms0[var]['mol'] = eta_i2/Tam*profs['Molecular_Backscatter_Coefficient'].profile*(pol[var][0]+pol[var][1]*pmol)
        ConstTerms0[var]['Ta'] = Tam
"""
# Example use of these constant terms
# #test_var = 'Raw_High_Gain_Total_Backscatter_Channel'
p_aer_test = 1.0-(profs['Particle_Depolarization'].profile)
#p_aer_test = 1.0
test_prof = {}
for test_var in raw_profs.keys():
    test_prof[test_var] = ConstTerms[test_var]['mult']*(profs['Aerosol_Backscatter_Coefficient'].profile*(pol[test_var][0]+p_aer_test*pol[test_var][1])+ConstTerms[test_var]['mol'])
"""

#%% Build TV 2D Image

range_offset = np.zeros(l_time.shape[0])

for piT in range(l_time.shape[0]-1):
    found_min = False
    cnt = 0
    dist = np.zeros(2)
    i_offset = 0    
    dist[0] = np.nanmean(np.sqrt((l_pos_x[0,:]-l_pos_x[piT+1,:])**2+(l_pos_y[0,:]-l_pos_y[piT+1,:])**2+(l_pos_z[0,:]-l_pos_z[piT+1,:])**2))
#    dist[0] = np.nanmean(np.sqrt((l_pos_x[piT,:]-l_pos_x[piT+1,:])**2+(l_pos_y[piT,:]-l_pos_y[piT+1,:])**2+(l_pos_z[piT,:]-l_pos_z[piT+1,:])**2))
    i_offset = 1
    dist[1] = np.nanmean(np.sqrt((l_pos_x[0,i_offset:]-l_pos_x[piT+1,:-i_offset])**2+(l_pos_y[0,i_offset:]-l_pos_y[piT+1,:-i_offset])**2+(l_pos_z[0,i_offset:]-l_pos_z[piT+1,:-i_offset])**2))
#    dist[1] = np.nanmean(np.sqrt((l_pos_x[piT,i_offset:]-l_pos_x[piT+1,:-i_offset])**2+(l_pos_y[piT,i_offset:]-l_pos_y[piT+1,:-i_offset])**2+(l_pos_z[piT,i_offset:]-l_pos_z[piT+1,:-i_offset])**2))

    if np.argmin(dist) == 0:
        # negative indexing
        i_offset = 1
        while(not found_min and i_offset < 1000):
            dist[1] = np.nanmean(np.sqrt((l_pos_x[0,:-i_offset]-l_pos_x[piT+1,i_offset:])**2+(l_pos_y[0,:-i_offset]-l_pos_y[piT+1,i_offset:])**2+(l_pos_z[0,:-i_offset]-l_pos_z[piT+1,i_offset:])**2))
#            dist[1] = np.nanmean(np.sqrt((l_pos_x[piT,:-i_offset]-l_pos_x[piT+1,i_offset:])**2+(l_pos_y[piT,:-i_offset]-l_pos_y[piT+1,i_offset:])**2+(l_pos_z[piT,:-i_offset]-l_pos_z[piT+1,i_offset:])**2))
            if dist[1] < dist[0]:
                dist[0] = dist[1]
            else:
                found_min = True
                range_offset[piT+1] = -(i_offset-1)
            i_offset+=1

    else:
        i_offset = 2
        while(not found_min and i_offset < 1000):
            dist[0] = np.nanmean(np.sqrt((l_pos_x[0,i_offset:]-l_pos_x[piT+1,:-i_offset])**2+(l_pos_y[0,i_offset:]-l_pos_y[piT+1,:-i_offset])**2+(l_pos_z[0,i_offset:]-l_pos_z[piT+1,:-i_offset])**2))
#            dist[0] = np.nanmean(np.sqrt((l_pos_x[piT,i_offset:]-l_pos_x[piT+1,:-i_offset])**2+(l_pos_y[piT,i_offset:]-l_pos_y[piT+1,:-i_offset])**2+(l_pos_z[piT,i_offset:]-l_pos_z[piT+1,:-i_offset])**2))
            if dist[0] < dist[1]:
                dist[1] = dist[0]
            else:
                found_min = True
                range_offset[piT+1] = i_offset-1
            i_offset+=1

#print(range_offset)
# I haven't validated this approach for desending aircraft
#range_offset = np.cumsum(range_offset)
range_offset = range_offset.astype(np.int)-np.int(range_offset.min())
range_offset = np.int(range_offset.max())-range_offset

             
range_array = np.zeros((l_pos_x.shape[0],np.int(l_pos_x.shape[1]+range_offset.max())))                
var = 'Raw_Molecular_Backscatter_Channel'

# define new profiles on the TV grid
rlen = l_pos_x.shape[1]
profs2D = {}
for var in profs.keys():
    profs2D[var]=profs[var].copy()
    profs2D[var].profile=np.zeros(range_array.shape)
    profs2D[var].profile_variance=np.zeros(range_array.shape)
    profs2D[var].range_array = range_array[0,:]  # dummy range array
    
raw_profs2D = {}
for var in raw_profs.keys():
    raw_profs2D[var]=raw_profs[var].copy()
    raw_profs2D[var].profile=np.zeros(range_array.shape)
    raw_profs2D[var].profile_variance=np.zeros(range_array.shape)
    raw_profs2D[var].range_array = range_array[0,:]  # dummy range array

ConstTerms = {}
for var in ConstTerms0.keys():
    ConstTerms[var] = {}
#    ConstTerms[var] = {'mult':np.zeros(range_array.shape),'mol':np.zeros(range_array.shape),'bg':ConstTerms[var]['bg'].copy(),'pol':ConstTerms[var]['pol'].copy(),'Ta':ConstTerms[var]['Ta']}
    for var2 in ConstTerms0[var].keys():
        if var2 == 'mult' or var2 == 'mol':
            ConstTerms[var][var2]=np.zeros(range_array.shape)
        else:
            try:
                ConstTerms[var][var2]=ConstTerms0[var][var2].copy()
            except AttributeError:
                ConstTerms[var][var2]=ConstTerms0[var][var2]
                
# create array containing ranges of the new rectangular grid
# and profiles of them
for piT in range(l_time.shape[0]-1):
    range_array[piT,range_offset[piT]:range_offset[piT]+rlen]=profs['Raw_Molecular_Backscatter_Channel'].range_array
    for var in profs.keys():
        profs2D[var].profile[piT,range_offset[piT]:range_offset[piT]+rlen]=profs[var].profile[piT,:]
        profs2D[var].profile_variance[piT,range_offset[piT]:range_offset[piT]+rlen]=profs[var].profile_variance[piT,:]
    for var in raw_profs.keys():
        raw_profs2D[var].profile[piT,range_offset[piT]:range_offset[piT]+rlen]=raw_profs[var].profile[piT,:]
        raw_profs2D[var].profile_variance[piT,range_offset[piT]:range_offset[piT]+rlen]=raw_profs[var].profile_variance[piT,:]
    for var in ConstTerms0.keys():
        ConstTerms[var]['mult'][piT,range_offset[piT]:range_offset[piT]+rlen]=ConstTerms0[var]['mult'][piT,:]
        ConstTerms[var]['mol'][piT,range_offset[piT]:range_offset[piT]+rlen]=ConstTerms0[var]['mol'][piT,:]

#%% Aerosol signal thresholding

"""
Attempt to remove noisy observations by duplicating the next lower signal with
a valid signal.
"""

## nan only thresholding
#AerBS0 = np.log(profs2D['Aerosol_Backscatter_Coefficient'].profile.copy())
#nan_count_last = 0
#nan_count = np.sum(np.isnan(AerBS0))
#while nan_count != nan_count_last:
#    inan = np.nonzero(np.isnan(AerBS0))
#    AerBS0[inan[0],inan[1]+1] = AerBS0[inan[0],inan[1]]
#    nan_count_last = nan_count
#    nan_count = np.sum(np.isnan(AerBS0))



## upper and lower thresholding
#Aer_threshold_lo = np.log(1e-9)  # backscatter coefficient threshold
#Aer_threshold_hi = np.log(1e-3)  # backscatter coefficient threshold
#
#AerBS0 = np.log(profs2D['Aerosol_Backscatter_Coefficient'].profile.copy())
#nan_count_last = 0
#nan_count = np.sum((AerBS0 < Aer_threshold_lo) + (AerBS0 > Aer_threshold_hi))
#while nan_count != nan_count_last:
#    inan = np.nonzero((AerBS0[:,1:] < Aer_threshold_lo)+(AerBS0[:,1:] >Aer_threshold_hi))
#    AerBS0[inan[0],inan[1]+1] = AerBS0[inan[0],inan[1]]
#    nan_count_last = nan_count
#    nan_count = np.sum((AerBS0[:,1:] < Aer_threshold_lo)+(AerBS0[:,1:] >Aer_threshold_hi))


# single threshold
Aer_threshold = np.log(1e-8)  # backscatter coefficient threshold

AerBS0 = np.log(profs2D['Aerosol_Backscatter_Coefficient'].profile.copy())
nan_count_last = 0
nan_count = np.sum(AerBS0 < Aer_threshold)
while nan_count != nan_count_last:
    inan = np.nonzero(AerBS0[:,1:] < Aer_threshold)
    AerBS0[inan[0],inan[1]+1] = AerBS0[inan[0],inan[1]]
    nan_count_last = nan_count
    nan_count = np.sum(AerBS0 < Aer_threshold)
    
AerBS0 = np.exp(AerBS0)

"""
# Show thresholding results
plt.figure()
plt.plot(AerBS0[5,:])
plt.plot(np.log(profs2D['Aerosol_Backscatter_Coefficient'].profile[5,:]))
"""

#MolCnts = profs2D['Raw_Molecular_Backscatter_Channel'].profile-ConstTerms['Raw_Molecular_Backscatter_Channel']['bg']
#ODest = np.cumsum(18*dR*profs2D['Aerosol_Backscatter_Coefficient'].profile,axis=1)
#Filt = np.cumsum((MolCnts < 24)*(ODest > 0.35),axis=1)>0



#for ai in range(Filt.shape[0]):
#    ifilt = np.nonzero(Filt[ai,:])[0]
#    if len(ifilt) > 1:
#        AerBS0[ai,ifilt] = AerBS0[ai,ifilt[0]]



#%% Denoising Script

"""
Denoising Main routine
"""
verify = True
verbose = True

#lam_array = np.array([10.0])
lam_array = np.array([np.nan])

dR = profs['Aerosol_Backscatter_Coefficient'].mean_dR

tdim = profs2D['Aerosol_Backscatter_Coefficient'].profile.shape[0]
rdim = profs2D['Aerosol_Backscatter_Coefficient'].profile.shape[1]

Interval = profs['Aerosol_Backscatter_Coefficient'].time.size   


fit_profs = {}
verify_profs = {}
for var in raw_profs.keys():
    # update to get actual photon counts in the profiles
    raw_profs2D[var].profile = raw_profs2D[var].profile*raw_profs[var].NumProfList[:,np.newaxis]
    raw_profs2D[var].profile_variance = raw_profs2D[var].profile_variance*raw_profs2D[var].NumProfList[:,np.newaxis]**2

    if verify:
        [fit_profs[var],verify_profs[var]] = raw_profs2D[var].p_thin()
        thin_adj = 0.5  # adjust gain on estimator if profiles have been poisson thinned
    else:
        fit_profs[var] = raw_profs2D[var].copy()
        thin_adj = 1.0
        
for var in raw_profs.keys():
    ConstTerms[var]['bg'] = ConstTerms[var]['bg']*thin_adj
      
#    
rate_adj = thin_adj*1.0/(fit_profs['Raw_Molecular_Backscatter_Channel'].shot_count*fit_profs['Raw_Molecular_Backscatter_Channel'].NumProfList*fit_profs['Raw_Molecular_Backscatter_Channel'].binwidth_ns*1e-9)


prefactor = 1.0# 1e-7

#LREstimate_buildConst_2D(geo_correct,MolScale,beta_m,range_array):
#ConstTerms = thin_adj*Nprof[:,np.newaxis]*LREstimate_buildConst_2D(geo_data,1,beta_mol.profile,MolRawE.range_array,MolRawE.mean_dt)
#ConstTerms = ConstTerms/beta_mol.profile



#xvalid2D = np.zeros(fit_profs['Raw_Molecular_Backscatter_Channel'].profile.shape)
#    CamList = np.zeros(Interval)
#    GmList = np.zeros(Interval)
#    GcList = np.zeros(Interval)
#    Dtmol = np.zeros(Interval)
#    Dtcomb = np.zeros(Interval)
#    fbgmol = np.zeros(Interval)
#    fbgcomb = np.zeros(Interval)
#sLR2D = np.zeros((tdim,rdim))
#beta_a_2D = np.zeros((tdim,rdim))
#alpha_a_2D = np.zeros((tdim,rdim))
#dep_a_2D = np.zeros((tdim,rdim))


#    fxVal = np.zeros(Interval)
opt_iterations = np.zeros(tdim)
opt_exit_mode = np.zeros(tdim)
#    str_exit_mode_list = []

# Find altitudes where there is molecular signal.  Only include those altitudes (but across all time) in the estimation
#xvalid = np.nonzero(np.sum(fit_profs['Raw_Molecular_Backscatter_Channel'].profile-fit_profs['Raw_Molecular_Backscatter_Channel'].bg[:,np.newaxis] \
#                        > 2 ,axis=0))[0]
#xvalid2D[:,xvalid] = 1    


errRecord = []
stepRecord = []

if verify:
    if np.isnan(lam_array).any():
        lam_array = np.logspace(-1,2,1)  # 47
else:
    lam_array = np.array([0])
#fitErrors = np.zeros(lam_array.size)
fitErrors = np.zeros(Num_reg_iter)
sol_List = []
lam_List = []
tv_list = []
#out_cond_array = np.zeros(lam_array.size)    
out_cond_array = np.zeros(Num_reg_iter)




lam_names = ['xB','xS','xP']  # variable names that have TV regularizer assigned to them

lam_sets = []
lam0 = {}
for lvar in lam_names:
    lam0[lvar] = 0

### Optimization Routine ###
#for i_lam in range(lam_array.size):  
for i_lam in range(Num_reg_iter):   
    
    if verify:
        lam = {}
        
        if Num_reg_iter > 1:
            for var in lam_range.keys():
                lam[var] = 10**(np.random.rand()*(lam_range[var][1]-lam_range[var][0])+lam_range[var][0])
        else:
            for var in lam_set.keys():
                lam[var] = lam_set[var]

        
#        lam['xB'] = 1e-3
#        lam['xS'] = 1e5
#        lam['xP'] = 1e5
        
#        lam['xB'] = 10**(np.random.rand()*1.0+1.5)
#        lam['xS'] = 10**(np.random.rand()*1.0-2.5)
#        lam['xP'] = 10**(np.random.rand()*1.0+0.0)
        
#        lam['xB'] = 69.9
#        lam['xS'] = 6.26e-3
#        lam['xP'] = 4.32
        
#        r1 = 10**(np.random.rand()*4.0-2.0)
#        r2 = 10**(np.random.rand()*4.0-2.0)
#        for lvar in lam_names:
#            if lvar == 'xS':
#                lam[lvar] = r2
#            else:
#                lam[lvar] = r1
                
#            lam[lvar] = 10**(np.random.rand()*4.0-2.0)
#            lam[lvar] = lam_array[i_lam]
#            lam[lvar] = 1e5

    else:
        lam = {}
        lam['xB'] = 10**1.0
        lam['xS'] = 10**1.5
        lam['xP'] = 10**1.0

    



    
    
#    bndsP = np.zeros((ix+2*xvalid.size*tdim,2))
#    bndsP[0,0] = 0.0
#    bndsP[0,1] = 0.1
#    bndsP[1,0] = 0.5
#    bndsP[1,1] = 1.1
#    bndsP[2,0] = 0.7
#    bndsP[2,1] = 1.5
#    
#    bndsP[3,0] = 0.0  # molecular deadtime
#    bndsP[3,1] = 100e-9
#    bndsP[4,0] = 0.0 # combined deeadtime
#    bndsP[4,1] = 100e-9
#    
#    bnds2D = np.zeros((tdim,2*xvalid.size))
#    bnds2D[:,:xvalid.size] = np.log(1.0*dR)  # lidar ratio lower limit
#    bnds2D[:,xvalid.size:] = np.log(1e-12)   # aerosol coefficient lower limit
#    bndsP[ix:,0] = bnds2D.flatten()
#    
#    bnds2D = np.zeros((tdim,2*xvalid.size))
#    bnds2D[:,:xvalid.size] = np.log(2e2*dR)  # lidar ratio upper limit
#    bnds2D[:,xvalid.size:] = np.log(1e-2)   # aerosol coefficient upper limit
#    bndsP[ix:,1] = bnds2D.flatten()
    
    # Gm, G hi/lo, Deadtimes on 4 channels
    x0 = {'xG':np.zeros(4),'xDT':np.zeros(4)}
    x0['xG'][0] = np.log(2.00)  # G mol
    x0['xG'][1] = np.log(1.0) # G hi
    x0['xG'][2] = np.log(0.004) # G lo
    x0['xG'][3] = np.log(1.0) # G cross
    x0['xDT'][0] = np.log(dead_time['Molecular_Backscatter_Channel'])  # molecular deadtime
    x0['xDT'][1] = np.log(dead_time['High_Gain_Total_Backscatter_Channel'])  # combined hi deadtime
    x0['xDT'][2] = np.log(dead_time['Low_Gain_Total_Backscatter_Channel'])  # combined lo deadtime
    x0['xDT'][3] = np.log(dead_time['Cross_Polarization_Channel'])  # cross deadtime
    
#    x02D[:,:rdim] = np.log(25.0-1)  # lidar ratio
#    x0['xS'] = np.log(np.random.rand(tdim,rdim)*30+20)  # lidar ratio
    x0['xS'] = 20*np.ones((tdim,rdim)) #  np.ones((tdim,rdim))*np.log(30.0-1)  # lidar ratio
    x0['xS'][np.nonzero(np.isnan(x0['xS']))] = opt_bounds['xS'][0] #np.log(1.0)  # get rid of invalid numbers

    #    x0['xB'] = np.log(profs2D['Aerosol_Backscatter_Coefficient'].profile)  # aerosol backscatter
    x0['xB'] = AerBS0.copy()  # aerosol backscatter
    x0['xB'][np.nonzero(np.isnan(x0['xB']))] = opt_bounds['xB'][0]   # get rid of invalid numbers
#    x0['xB'][np.nonzero(x0['xB']<1e-9)] = 1e-9    # get rid of invalid numbers
#    x0['xB'][np.nonzero(x0['xB']>1e-2)] = 1e-2    # get rid of invalid numbers
    x0['xP'] = 1-profs2D['Particle_Depolarization'].profile #np.tan((0.5-profs2D['Particle_Depolarization'].profile)*np.pi)  # aerosol depolarization
    x0['xP'][np.nonzero(np.isnan(x0['xP']))] = opt_bounds['xP'][1] #np.tan((0.5-0.05)*np.pi)    # get rid of invalid numbers
#    x0['xP'][np.nonzero(x0['xP'] > 0.99)] = 0.99 #np.tan((0.5-0.01)*np.pi)    # get rid of invalid numbers
#    x0['xP'][np.nonzero(x0['xP'] < 0.2)] = 0.2 #np.tan((0.5-0.01)*np.pi)    # get rid of invalid numbers
#    x0['xP'][np.nonzero(x0['xP'] >= 32)] = np.tan((0.5-0.01)*np.pi)    # get rid of invalid numbers
    
    
#    condition_functions = {'xB':lambda x,y: mle.cond_exp(x,np.nanstd(np.log(x0['xB'])),-np.nanmean(np.log(x0['xB'])),operation=y),
#                       'xS':lambda x,y: mle.cond_linear(x,40,-20,operation=y),
#                       'xP':lambda x,y: mle.cond_arctan(x,np.nanstd(x0['xP']),-np.nanmean(x0['xP']),0,1.0,operation=y)}    
    
    condition_functions = {'xB':lambda x,y: mle.cond_exp(x,np.log(opt_bounds['xB'][1])-np.log(opt_bounds['xB'][0]),np.log(opt_bounds['xB'][0])/(np.log(opt_bounds['xB'][1])-np.log(opt_bounds['xB'][0])),operation=y),
                       'xS':lambda x,y: mle.cond_linear(x,opt_bounds['xS'][1]-opt_bounds['xS'][0],1/(opt_bounds['xS'][1]-opt_bounds['xS'][0]),operation=y),
                       'xP':lambda x,y: mle.cond_linear(x,1.0,0.0,operation=y)}   
    
#    bnds = {'xS':sorted([condition_functions['xS'](1,'inverse'),condition_functions['xS'](100,'inverse')]),
#            'xB':sorted([condition_functions['xB'](1e-9,'inverse'),condition_functions['xB'](1e-2,'inverse')]),
#            'xP':sorted([condition_functions['xP'](0,'inverse'),condition_functions['xP'](1,'inverse')])} 
    
    # apply condition functions to initial conditions and variable bounds
    bnds = {}
    for var in condition_functions:
        x0[var] = condition_functions[var](x0[var],'inverse')
        bnds[var] = sorted([condition_functions[var](opt_bounds[var][0],'inverse'),condition_functions[var](opt_bounds[var][1],'inverse')])
        x0[var][x0[var] < bnds[var][0]] = bnds[var][0]
        x0[var][x0[var] > bnds[var][1]] = bnds[var][1]
        
    FitProf = lambda x: mle.GVHSRL_sparsa_Error(x,fit_profs,ConstTerms,lam,dt=rate_adj,cond_fun=condition_functions)

    FitProfDeriv = lambda x: mle.GVHSRL_sparsa_Error_Gradient(x,fit_profs,ConstTerms,lam,dt=rate_adj,cond_fun=condition_functions)
    
#    scale = {}
#    scale['xS'] = 40 #np.abs(np.log(40.0-1)) #np.nanmean(np.abs(x0['xS']))
#    scale['xB'] = np.nanmean(np.abs(x0['xB']))  # 7.0
#    scale['xP'] = 10.0 #350; #np.nanmean(np.abs(x0['xP']))
    
#    bnds = {'xS':sorted([1.0/scale['xS'],100.0/scale['xS']]),
#            'xB':sorted([np.log(1e-9/scale['xB']),np.log(1e-2/scale['xB'])]),
#            'xP':sorted([0,1.0/scale['xP']])}    
    
#    for var in scale.keys():
#        x0[var] = x0[var]/scale[var]
        
#    FitProf = lambda x: mle.GVHSRL_sparsa_Error(x,fit_profs,ConstTerms,lam,dt=rate_adj,scale=scale)
#
#    FitProfDeriv = lambda x: mle.GVHSRL_sparsa_Error_Gradient(x,fit_profs,ConstTerms,lam,dt=rate_adj,scale=scale)
    
    
#    prof_sol = Build_GVHSRL_Profiles(x0,ConstTerms,dt=rate_adj,ix=ix,return_params=True)
    
    """
    #Test code for Gradient
    gradAnalytic = FitProfDeriv(x0)
    gradNum = mle.Num_Gradient_Dict(FitProf,x0)
    for var in gradNum.keys():
        plt.figure()
        plt.plot(gradNum[var].flatten())
        plt.plot(gradAnalytic[var].flatten(),'--')
        plt.title(str(var))
    """
    
    
    """
    Test Code:
    var = 'Raw_Molecular_Backscatter_Channel'; 
    plt.figure(); 
    plt.semilogy(prof_sol[var][40,:]); 
    plt.semilogy(fit_profs[var].profile[40,:])
    
    FitProf(x0)
    """    
    
 
    
    sol,[error_hist,step_hist]= mle.GVHSRL_sparsa_optimizor(FitProf,FitProfDeriv,x0,lam,sub_eps=1e-5,step_eps=step_eps,opt_cnt_min=10,opt_cnt_max=max_iter,cnt_alpha_max=10,sigma=1e-5,verbose=False,alpha = 1e5,bnds=bnds)
    
    if not verify:
        prof_sol = mle.Build_GVHSRL_sparsa_Profiles(sol,ConstTerms,dt=rate_adj,return_params=True,cond_fun=condition_functions)

    
    ProfileLogError = mle.GVHSRL_sparsa_Error(sol,verify_profs,ConstTerms,lam0,dt=rate_adj)    
    

    errRecord.extend([error_hist])
    stepRecord.extend([step_hist])
    lam_List.extend([lam.copy()])
    lam_sets.extend([[lam['xB'],lam['xS'],lam['xP'],ProfileLogError]])
    tv_sublist = []
    for xvar in ['xB','xS','xP']:
        tv_sublist.extend([np.nansum(np.diff(sol[xvar],axis=1))+np.nansum(np.diff(sol[xvar],axis=0))])
    tv_list.extend([tv_sublist])
    fitErrors[i_lam] = np.nansum(ProfileLogError)
    sol_List.extend([sol.copy()])        
    if verbose:
        print('Iteration: %d'%i_lam)
        print('Log Error: %f'%fitErrors[i_lam])
        print('lambda B, S, P')
        print('    %e | %e |  %e'%(lam['xB'],lam['xS'],lam['xP']))
        print('')

### End Optimization Routine ###

lam_sets = np.array(lam_sets)
tv_list = np.array(tv_list)
xG = sol['xG']
xDT = sol['xDT']
xB = sol['xB']
xS = sol['xS']
xP = sol['xP']
if save_results:
    np.savez(savefigpath+'_3D_opt_results_'+datetime.datetime.now().strftime('%Y%m%dT%H%M'),lam_sets=lam_sets,tv_list=tv_list,xG=xG,xDT=xDT,xB=xB,xS=xS,xP=xP,step_eps=step_eps)

## 1D regularizer
#plt.figure()
#plt.semilogx(lam_array,fitErrors)
#plt.xlabel('Regularizer')
#plt.ylabel('Fit Error')
#plt.grid(b=True)
#if save_results:
#    plt.savefig(savefigpath+'Regularizer_GV_HSRL.png',dpi=300)

## 2D regularizer
#plt.figure()
#plt.scatter(lam_sets[:,0],lam_sets[:,1],c=lam_sets[:,3])
#plt.xlabel(r'$\lambda_{\beta}$')
#plt.ylabel(r'$\lambda_{s}$')
#plt.xscale('log')
#plt.yscale('log')
#plt.colorbar()

# 3D regularizer
fig = plt.figure()
mx = fig.add_subplot(111,projection='3d')
mx.set_xlabel(r'$\log_{10} \lambda_{\beta}$')
mx.set_ylabel(r'$\log_{10} \lambda_{s}$')
mx.set_zlabel(r'$\log_{10} \lambda_{p}$')
pdata = mx.scatter(np.log10(lam_sets[:,0]),np.log10(lam_sets[:,1]),np.log10(lam_sets[:,2]),c=lam_sets[:,3])
if save_results:
    plt.savefig(savefigpath+'_RegularizerSpace_GV_HSRL.png',dpi=300)

isol = np.argmin(fitErrors)

print('Solution:')
print('lambda B, S, P')
print('    %e | %e |  %e'%(lam_sets[isol,0],lam_sets[isol,1],lam_sets[isol,2]))

plt.figure()
plt.plot(errRecord[isol])
plt.xlabel('Iterations')
plt.ylabel('Log-likelihood')
plt.grid(b=True)
if save_results:
    plt.savefig(savefigpath+'_Loglikelihood_GV_HSRL.png',dpi=300)

plt.figure()
plt.semilogy(stepRecord[isol][1:])
plt.xlabel('Iterations')
plt.ylabel('Step Evaluation')
plt.grid(b=True)
if save_results:
    plt.savefig(savefigpath+'_StepEvaluation_GV_HSRL.png',dpi=300)
    
plt.figure(); 
plt.plot([len(x) for x in errRecord],[x[-1] for x in errRecord],'.')
plt.xlabel('Number of Iterations')
plt.ylabel('LogLikelihood')
plt.grid(b=True)

prof_sol = mle.Build_GVHSRL_sparsa_Profiles(sol_List[isol],ConstTerms,dt=rate_adj,return_params=True,cond_fun=condition_functions)

prof_x0 = mle.Build_GVHSRL_sparsa_Profiles(x0,ConstTerms,dt=rate_adj,return_params=True,cond_fun=condition_functions)

tind = np.int(np.round(np.random.rand()*sol['xB'].shape[0]))-1
subnum = 410
plt.figure()
for ifig,var in enumerate(fit_profs.keys()):
    plt.subplot(subnum+1+ifig)
    plt.semilogy(fit_profs[var].profile[tind,:])
    plt.semilogy(verify_profs[var].profile[tind,:])
    plt.semilogy(prof_sol[var][tind,:])
if save_results:
    plt.savefig(savefigpath+'Fit_%d_index_GV_HSRL.png'%tind,dpi=300)

plt.figure()
plt.subplot(141)
range_mask = np.ones(profs2D['Aerosol_Backscatter_Coefficient'].range_array.size)
range_mask[profs2D['Aerosol_Backscatter_Coefficient'].range_array==0] = np.nan
plt.semilogx(profs2D['Aerosol_Backscatter_Coefficient'].profile[tind,:],profs2D['Aerosol_Backscatter_Coefficient'].range_array*range_mask)
plt.semilogx(prof_sol['Backscatter_Coefficient'][tind,:],profs2D['Aerosol_Backscatter_Coefficient'].range_array*range_mask)
plt.xlim([1e-7,1e-2])
plt.xlabel(r'$\beta_a$')
plt.ylabel(r'Range [m]')
plt.grid(b=True)
plt.subplot(142)
plt.semilogx(profs2D['Aerosol_Extinction_Coefficient'].profile[tind,:],profs2D['Aerosol_Extinction_Coefficient'].range_array*range_mask)
plt.semilogx(prof_sol['Backscatter_Coefficient'][tind,:]*prof_sol['Lidar_Ratio'][tind,:],profs2D['Aerosol_Extinction_Coefficient'].range_array*range_mask)
plt.xlim([1e-5,1e0])
plt.xlabel(r'$\alpha_a$')
plt.grid(b=True)
plt.subplot(143)
plt.plot(profs2D['Aerosol_Extinction_Coefficient'].profile[tind,:]/profs2D['Aerosol_Backscatter_Coefficient'].profile[tind,:],profs2D['Aerosol_Extinction_Coefficient'].range_array*range_mask)
plt.plot(prof_sol['Lidar_Ratio'][tind,:],profs2D['Aerosol_Extinction_Coefficient'].range_array*range_mask)
plt.xlim([0,50])
plt.xlabel(r'$s_{LR}$')
plt.grid(b=True)
plt.subplot(144)
plt.plot(profs2D['Particle_Depolarization'].profile[tind,:],profs2D['Particle_Depolarization'].range_array*range_mask)
plt.plot(1-prof_sol['Polarization'][tind,:],profs2D['Particle_Depolarization'].range_array*range_mask)
plt.xlim([0,1])
plt.xlabel(r'$d_{a}$')
plt.grid(b=True)

if save_results:
    plt.savefig(savefigpath+'_FitParams_%d_index_GV_HSRL.png'%tind,dpi=300)

#plt.show()


denoise_labels = ['Aerosol_Backscatter_Coefficient','Particle_Depolarization','Aerosol_Extinction_Coefficient','Lidar_Ratio']
denoise_profs = {}
for var in denoise_labels:
    try:    
        denoise_profs[var] = profs[var].copy()
        denoise_profs[var].label = 'Denoised '+denoise_profs[var].label
        denoise_profs[var].descript = 'PTV Denoised ' + denoise_profs[var].descript
    except KeyError:
        denoise_profs[var] = profs[denoise_labels[0]].copy()
        denoise_profs[var].label = 'Denoised Lidar Ratio'
        denoise_profs[var].profile_type = 'sr'
        denoise_profs[var].descript = 'PTV Denoised Aerosol Lidar Ratio'

# create array containing ranges of the new rectangular grid
# and profiles of them
for piT in range(l_time.shape[0]-1):
    for var in denoise_profs.keys():
        if 'Backscatter' in var:
            denoise_profs[var].profile[piT,:] = prof_sol['Backscatter_Coefficient'][piT,range_offset[piT]:range_offset[piT]+rlen]
        elif 'Extinction' in var:
            denoise_profs[var].profile[piT,:] = prof_sol['Backscatter_Coefficient'][piT,range_offset[piT]:range_offset[piT]+rlen]*prof_sol['Lidar_Ratio'][piT,range_offset[piT]:range_offset[piT]+rlen]
        elif 'Lidar_Ratio' in var:
            denoise_profs[var].profile[piT,:] = prof_sol['Lidar_Ratio'][piT,range_offset[piT]:range_offset[piT]+rlen]
        elif 'Depolarization' in var:
            denoise_profs[var].profile[piT,:] = 1-prof_sol['Polarization'][piT,range_offset[piT]:range_offset[piT]+rlen]
            

for plt_var in denoise_profs.keys():
    lplt.scatter_z(denoise_profs[plt_var],lidar_pointing=lidar_data['lidar_pointing'],
               lidar_alt=aircraft_data['GGALT'],climits=plot_settings[plt_var]['climits'],
               cmap=plot_settings[plt_var]['colormap'],scale=plot_settings[plt_var]['scale'],s=2,t_axis_scale=5.0)
    if save_results:
        plt.savefig(savefigpath+plt_var+'_PTV_GV_HSRL.png',dpi=600)
        
    if plt_var in profs.keys():
        lplt.scatter_z(profs[plt_var],lidar_pointing=lidar_data['lidar_pointing'],
               lidar_alt=aircraft_data['GGALT'],climits=plot_settings[plt_var]['climits'],
               cmap=plot_settings[plt_var]['colormap'],scale=plot_settings[plt_var]['scale'],s=2,t_axis_scale=5.0)
        
        if save_results:
            plt.savefig(savefigpath+'_'+plt_var+'_GV_HSRL.png',dpi=600)

plt.show()

if save_results:
    for var in denoise_profs.keys():
        denoise_profs[var].write2nc(save_path+savencfile)
    
    for var in profs.keys():
        profs[var].write2nc(save_path+savencfile)
    
    for var in lidar_data.keys():
        lp.write_var2nc(lidar_data[var],var,save_path+savencfile)
        
    for var in aircraft_data.keys():
        lp.write_var2nc(aircraft_data[var],var,save_path+savencfile)
        
    lp.write_var2nc(lam_sets,'regularizer_data',save_path+savencfile,description='Regularizer values organized by [Backscatter, Lidar Ratio, Depolarization, Log-Likelihood]')
    lp.write_var2nc(tv_list,'tv_data',save_path+savencfile,description='Profile total variation organized by [Backscatter, Lidar Ratio, Depolarization]')


