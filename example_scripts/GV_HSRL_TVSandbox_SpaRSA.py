# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 08:20:36 2018

@author: mhayman
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

import LidarProfileFunctions as lp
import LidarPlotFunctions as lplt
import GVHSRLlib as gv
import MLELidarProfileFunctions as mle

import datetime

import netCDF4 as nc4

import json


    
    
    
# path to files being analyzed
#ncpath = '/scr/sci/mhayman/SOCRATES/range_centered/test_data/SOCRATESrf03/'
ncpath = '/Users/mhayman/Documents/HSRL/GVHSRL/test_data/SOCRATESrf03/'
# list of files to analyze
nclist0 = ['SOCRATESrf03_GVHSRL_20180123T0210_20180123T0220.nc']

# if desired, loaded files will be constrainted by time
time_start = datetime.datetime(year=1900,month=1,day=1)
time_stop = datetime.datetime(year=2100,month=1,day=1)

load_mask = False

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
        'climits':[0,0.6],
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
        'climits':[10,50],
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
#%%  Plot Sample Profile

var = 'Aerosol_Backscatter_Coefficient'
lplt.scatter_z(profs[var],lidar_pointing=lidar_data['lidar_pointing'],
               lidar_alt=aircraft_data['GGALT'],climits=plot_settings[var]['climits'],
               cmap=plot_settings[var]['colormap'],scale=plot_settings[var]['scale'],s=2,t_axis_scale=5.0)

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

t_lim = [26*3600+13*60,26*3600+13.4*60]
r_lim = [0,2e3]

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

#cal_file_path = '/h/eol/mhayman/PythonScripts/HSRL_Processing/GV-HSRL-Python/calibrations/cal_files/'
cal_file_path = '/Users/mhayman/Documents/Python/Lidar/GV-HSRL-Python/calibrations/cal_files/'
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

#%% Denoising Script

"""
Denoising Main routine
"""
verify = True
verbose = False

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

if verify:
    if np.isnan(lam_array).any():
        lam_array = np.logspace(-3,3,100)  # 47
else:
    lam_array = np.array([0])
fitErrors = np.zeros(lam_array.size)
sol_List = []
lam_List = []
out_cond_array = np.zeros(lam_array.size)    

lam_names = ['xB','xS','xP']  # variable names that have TV regularizer assigned to them

lam_sets = []
lam0 = {}
for lvar in lam_names:
    lam0[lvar] = 0

### Optimization Routine ###
for i_lam in range(lam_array.size):    
    
    if verify:
        lam = {}
        r1 = 10**(np.random.rand()*1)
        r2 = 10**(np.random.rand()*0.5+2)
        for lvar in lam_names:
            if lvar == 'xS':
                lam[lvar] = r2
            else:
                lam[lvar] = r1
                
#            lam[lvar] = 10**(np.random.rand()*6-3)
#            lam[lvar] = lam_array[i_lam]

    else:
        lam = {}
        lam['xB'] = 10**1.0
        lam['xS'] = 10**1.5
        lam['xP'] = 10**1.0

    
    #nonlinear
    # GVHSRL_FitError(x,fit_profs,Const,lam,meshTV,dt=1.0,ix=8,weights=np.array([1])):
    FitProf = lambda x: mle.GVHSRL_sparsa_Error(x,fit_profs,ConstTerms,lam,dt=rate_adj)
    #GVHSRL_FitError(x,fit_profs,ConstTerms,lam,mesh_pts,dt=rate_adj)
    FitProfDeriv = lambda x: mle.GVHSRL_sparsa_Error_Gradient(x,fit_profs,ConstTerms,lam,dt=rate_adj)
    #GVHSRL_FitError_Gradient(x,fit_profs,ConstTerms,lam,mesh_pts,dt=rate_adj)
    
    
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
    x0['xS'] = np.ones((tdim,rdim))*np.log(25.0-1)  # lidar ratio
    x0['xS'][np.nonzero(np.isnan(x0['xS']))] = np.log(1.0)  # get rid of invalid numbers

    x0['xB'] = np.log(profs2D['Aerosol_Backscatter_Coefficient'].profile)  # aerosol backscatter
    x0['xB'][np.nonzero(np.isnan(x0['xB']))] = np.log(1e-12)    # get rid of invalid numbers
    x0['xB'][np.nonzero(x0['xB']<np.log(1e-12))] = np.log(1e-12)    # get rid of invalid numbers
    x0['xP'] = np.tan((0.5-profs2D['Particle_Depolarization'].profile)*np.pi)  # aerosol depolarization
    x0['xP'][np.nonzero(np.isnan(x0['xP']))] = -1e4    # get rid of invalid numbers
    
    
#    prof_sol = Build_GVHSRL_Profiles(x0,ConstTerms,dt=rate_adj,ix=ix,return_params=True)
    
    """
    #Test code for Gradient
    gradNum = mle.Num_Gradient(FitProf,x0)
    plt.figure()
    plt.plot(gradNum)
    plt.plot(FitProfDeriv(x0),'--')
    """
    
    
    """
    Test Code:
    var = 'Raw_Molecular_Backscatter_Channel'; 
    plt.figure(); 
    plt.semilogy(prof_sol[var][40,:]); 
    plt.semilogy(fit_profs[var].profile[40,:])
    
    FitProf(x0)
    """    
    
 
    
    sol,error_hist= mle.GVHSRL_sparsa_optimizor(FitProf,FitProfDeriv,x0,lam,sub_eps=1e-5,step_eps=1e-30,opt_cnt_max=100,cnt_alpha_max=10,sigma=1e-5,verbose=False,alpha = 1e20)
    
    if not verify:
        prof_sol = mle.Build_GVHSRL_sparsa_Profiles(sol,ConstTerms,dt=rate_adj,return_params=True)

    
    ProfileLogError = mle.GVHSRL_sparsa_Error(sol,verify_profs,ConstTerms,lam0,dt=rate_adj)    
    

    errRecord.extend([error_hist])
    lam_List.extend([lam.copy()])
    lam_sets.extend([[lam['xB'],lam['xS'],lam['xP'],ProfileLogError]])
    fitErrors[i_lam] = np.nansum(ProfileLogError)
    sol_List.extend([sol.copy()])        
    if verbose:
        print('Log Error: %f'%fitErrors[i_lam])

### End Optimization Routine ###

#plt.figure()
#plt.semilogx(lam_array,fitErrors)
#plt.xlabel('Regularizer')
#plt.ylabel('Fit Error')
#plt.grid(b=True)

lam_sets = np.array(lam_sets)
plt.figure()
plt.scatter(lam_sets[:,0],lam_sets[:,1],c=lam_sets[:,3])
plt.xlabel(r'$\lambda_{\beta}$')
plt.ylabel(r'$\lambda_{s}$')
plt.xscale('log')
plt.yscale('log')
plt.colorbar()

isol = np.argmin(fitErrors)

plt.figure()
plt.plot(errRecord[isol])
plt.xlabel('Iterations')
plt.ylabel('Log-likelihood')
plt.grid(b=True)

prof_sol = mle.Build_GVHSRL_sparsa_Profiles(sol_List[isol],ConstTerms,dt=rate_adj,return_params=True)

tind = np.int(np.round(np.random.rand()*sol['xB'].shape[0]))
subnum = 410
plt.figure()
for ifig,var in enumerate(fit_profs.keys()):
    plt.subplot(subnum+1+ifig)
    plt.semilogy(fit_profs[var].profile[tind,:])
    plt.semilogy(verify_profs[var].profile[tind,:])
    plt.semilogy(prof_sol[var][tind,:])



plt.figure()
plt.subplot(411)
plt.semilogy(profs2D['Aerosol_Backscatter_Coefficient'].profile[tind,:])
plt.semilogy(prof_sol['Backscatter_Coefficient'][tind,:])
plt.ylim([1e-9,1e-3])
plt.subplot(412)
plt.semilogy(profs2D['Aerosol_Extinction_Coefficient'].profile[tind,:])
plt.semilogy(prof_sol['Backscatter_Coefficient'][tind,:]*prof_sol['Lidar_Ratio'][tind,:])
plt.ylim([1e-7,1e-1])
plt.subplot(413)
plt.plot(profs2D['Aerosol_Extinction_Coefficient'].profile[tind,:]/profs2D['Aerosol_Backscatter_Coefficient'].profile[tind,:])
plt.plot(prof_sol['Lidar_Ratio'][tind,:])
plt.ylim([0,100])
plt.subplot(414)
plt.plot(profs2D['Particle_Depolarization'].profile[tind,:])
plt.plot(1-prof_sol['Polarization'][tind,:])
plt.ylim([0,1])

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
        elif 'Lidar Ratio' in var:
            denoise_profs[var].profile[piT,:] = prof_sol['Lidar_Ratio'][piT,range_offset[piT]:range_offset[piT]+rlen]
        elif 'Depolarization' in var:
            denoise_profs[var].profile[piT,:] = 1-prof_sol['Polarization'][piT,range_offset[piT]:range_offset[piT]+rlen]
            

for plt_var in denoise_profs.keys():
    lplt.scatter_z(denoise_profs[plt_var],lidar_pointing=lidar_data['lidar_pointing'],
               lidar_alt=aircraft_data['GGALT'],climits=plot_settings[plt_var]['climits'],
               cmap=plot_settings[plt_var]['colormap'],scale=plot_settings[plt_var]['scale'],s=2,t_axis_scale=5.0)
    if plt_var in profs.keys():
        lplt.scatter_z(profs[plt_var],lidar_pointing=lidar_data['lidar_pointing'],
               lidar_alt=aircraft_data['GGALT'],climits=plot_settings[plt_var]['climits'],
               cmap=plot_settings[plt_var]['colormap'],scale=plot_settings[plt_var]['scale'],s=2,t_axis_scale=5.0)

plt.show()



#beta_a_denoise = profs['Aerosol_Backscatter_Coefficient'].copy()
#beta_a_denoise.profile = prof_sol['Backscatter_Coefficient'].copy()
#beta_a_denoise.label = 'Denoised ' + beta_a_denoise.label
#
#s_LR_denoise = profs['Aerosol_Backscatter_Coefficient'].copy()
#s_LR_denoise.profile = prof_sol['Lidar_Ratio'].copy()
#s_LR_denoise.label = 'Denoised Lidar Ratio'
#s_LR_denoise.profile_type = '$sr$'
#
#dPart_denoise = profs['Particle_Depolarization'].copy()
#dPart_denoise.profile = 2*(1-prof_sol['Polarization'].copy())
#dPart_denoise.label = 'Denoised ' + dPart_denoise.label
#
#
#var = 'Aerosol_Backscatter_Coefficient'
#lplt.scatter_z(profs[var],lidar_pointing=lidar_data['lidar_pointing'],
#               lidar_alt=aircraft_data['GGALT'],climits=plot_settings[var]['climits'],
#               cmap=plot_settings[var]['colormap'],scale=plot_settings[var]['scale'],s=2,t_axis_scale=5.0)
#
#lplt.scatter_z(beta_a_denoise,lidar_pointing=lidar_data['lidar_pointing'],
#               lidar_alt=aircraft_data['GGALT'],climits=plot_settings[var]['climits'],
#               cmap=plot_settings[var]['colormap'],scale=plot_settings[var]['scale'],s=2,t_axis_scale=5.0)
#
#var = 'Lidar_Ratio'               
#lplt.scatter_z(s_LR_denoise,lidar_pointing=lidar_data['lidar_pointing'],
#               lidar_alt=aircraft_data['GGALT'],climits=plot_settings[var]['climits'],
#               cmap=plot_settings[var]['colormap'],scale=plot_settings[var]['scale'],s=2,t_axis_scale=5.0)
#
#var = 'Particle_Depolarization'
#lplt.scatter_z(dPart_denoise,lidar_pointing=lidar_data['lidar_pointing'],
#               lidar_alt=aircraft_data['GGALT'],climits=plot_settings[var]['climits'],
#               cmap=plot_settings[var]['colormap'],scale=plot_settings[var]['scale'],s=2,t_axis_scale=5.0)
 
"""
    
fit_mol_2D,fit_comb_2D = ProfilesTotalxvalid_2D(sol_List[isol],xvalid,tdim,beta_mol.profile,ConstTerms,Mprof_bg=FitMol_bg,Cprof_bg=FitComb_bg,dt=rate_adj)

Cam = sol_List[isol][0]
Gm = sol_List[isol][1]
Gc = sol_List[isol][2]      
if print_sol and verify:
    print('Final Solution:')
    print('lam: %f\nCam = %f\nGm = %f\nGc = %f\nMol DT = %f ns\nComb DT = %f ns\nFitError:%f'%(lam_array[isol],Cam,Gm,Gc,sol_List[isol][3]*1e9,sol_List[isol][4]*1e9,fitErrors[isol]))
    print('Output Flag: %d'%out_cond_array[isol])
    print('Output Flag Definition: %s'%scipy.optimize.tnc.RCSTRINGS[out_cond_array[isol]])

if plotfigs:
    if verify:
        plt.figure()
        plt.semilogx(lam_array,fitErrors)
        plt.ylabel('Log Error')
        plt.xlabel('Fit Index')       

sol2D = sol_List[isol][ix:].reshape((tdim,2*xvalid.size))
beta_a_2D[:,xvalid] = np.exp(sol2D[:,xvalid.size:])
sLR2D[:,xvalid] = np.exp(sol2D[:,:xvalid.size])/dR

#    Cam = sol1D[0]
#    Gm = sol1D[1]
#    Gc = sol1D[2]
#    DTmol = sol1D[3]
#    DTcomb = sol1D[4]

cal_params = {'Cam':sol_List[isol][0],'Gm':sol_List[isol][1],'Gc':sol_List[isol][2],'DTmol':sol_List[isol][3],'DTcomb':sol_List[isol][4]}


#    fit_mol_2D[pI,:],fit_comb_2D[pI,:] = ProfilesTotalxvalid(wMol,xvalid,beta_m_sonde.profile[pI,:],ConstTerms[pI,:],Mprof_bg=FitMol_bg[pI],Cprof_bg=FitComb_bg[pI])
#ProfilesTotalxvalid_2D(x,xvalid,tdim,mol_bs_coeff,Const,Mprof_bg=0,Cprof_bg=0,dt=np.array([1.0]))


alpha_a_2D = beta_a_2D*sLR2D



#    # attempt to remove profiles where the fit was bad
#    pbad = np.nonzero(np.abs(np.diff(ProfileError)/ProfileError[:-1])>3)[0]+1
#    xvalid2D[pbad,:]= 0
#    pbad = np.nonzero(np.abs(np.diff(ProfileError)/ProfileError[1:])>3)[0]
#    xvalid2D[pbad,:]= 0


# Merge with original aerosol data set
beta_merge = beta_aer.copy()
beta_a_2D_adj = beta_a_2D[:,:beta_aer.profile.shape[1]]
iMLE = np.nonzero(xvalid2D)
beta_merge.profile[iMLE] = beta_a_2D_adj[iMLE]
beta_merge.descript = 'Maximum Likelihood Estimate of Aerosol Backscatter Coefficient in m^-1 sr^-1'

sLR_mle = beta_aer.copy()
sLR_mle.profile = sLR2D[:,:beta_aer.profile.shape[1]].copy()
sLR_mle.descript = 'Maximum Likelihood Estimate of Aerosol Lidar Ratio sr'
sLR_mle.label = 'Aerosol Lidar Ratio'
sLR_mle.profile_type = '$sr$'
sLR_mle.profile_variance = sLR_mle.profile_variance*0.0

beta_a_mle = beta_aer.copy()
beta_a_mle.profile = beta_a_2D[:,:beta_aer.profile.shape[1]].copy()
beta_a_mle.descript = 'Maximum Likelihood Estimate of Aerosol Backscatter Coefficient in m^-1 sr^-1'
beta_a_mle.label = 'Aerosol Backscatter Coefficient'
beta_a_mle.profile_type = '$m^{-1}sr^{-1}$'
beta_a_mle.profile_variance = beta_a_mle.profile_variance*0.0

alpha_a_mle = beta_aer.copy()
alpha_a_mle.profile = alpha_a_2D[:,:beta_aer.profile.shape[1]].copy()
alpha_a_mle.descript = 'Maximum Likelihood Estimate of Aerosol Extinction Coefficient in m^-1'
alpha_a_mle.label = 'Aerosol Extinction Coefficient'
alpha_a_mle.profile_type = '$m^{-1}$'
alpha_a_mle.profile_variance = alpha_a_mle.profile_variance*0.0

fit_mol_mle = MolRawE.copy()
fit_mol_mle.profile = fit_mol_2D.copy()
fit_mol_mle.slice_range_index(range_lim=[0,beta_aer.profile.shape[1]])
fit_mol_mle.descript = 'Maximum Likelihood Estimate of ' + fit_mol_mle.descript
fit_comb_mle = CombRawE.copy()
fit_comb_mle.profile = fit_comb_2D.copy()
fit_comb_mle.slice_range_index(range_lim=[0,beta_aer.profile.shape[1]])
fit_comb_mle.descript = 'Maximum Likelihood Estimate of ' + fit_comb_mle.descript

xvalid_mle = beta_aer.copy()
xvalid_mle.profile = xvalid2D[:,:beta_aer.profile.shape[1]].copy()
xvalid_mle.descript = 'Maximum Likelihood Estimated Data Points'
xvalid_mle.label = 'MLE Data Point Mask'
xvalid_mle.profile_type = 'MLE Data Point Mask'
xvalid_mle.profile_variance = xvalid_mle.profile*0


#    if use_mask:
#        NanMask = np.logical_or(Molecular.profile < 4.0,CombHi.profile < 4.0)
#        beta_merge.profile = np.ma.array(beta_merge,mask=NanMask)


#return beta_merge,beta_a_mle,sLR_mle,alpha_a_mle,fit_mol_mle,fit_comb_mle,xvalid_mle,ProfileLogError,cal_params

"""