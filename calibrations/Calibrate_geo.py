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
    airborne_in = input('Is this an airborne calibration? [y/n]')
    if airborne_in == 'y' or airborne_in == 'Y':
#        start_hr = np.float(input("Start Hour (UTC): "))
#        stop_hr = np.float(input("Duration (hours): "))
        airborne = True
    else:
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
    
    'airborne':False, # is the lidar airborne (downward poinging) 
    'Airspeed_Threshold':15, # threshold for determining start and end of the flight (in m/s)
    'ground_range_buffer':3e3,  # max altitude of ground [m] to avoid counting it as a cloud    
    
    'pol_xtalk':0.015,
    
    
    'FilterI2':False,  # only include data where the I2 cell is removed
    
    'diff_geo_correct':True,
    'hsrl_rb_adjust':True,
    'RemoveCals':True, # removes instances where the I2 cell is removed
    'load_reanalysis':True,
    'use_bs_thr':True,
    
    'EstimateExtinction':True,  # attempt to apply an extinction correction
    'ExtinctionAlt':5e3,     # max altitude to apply extinction correction
    
    'sg_win':11,   # window size in SG filter
    'sg_order':5  # order of SG filter

    
    }

if airborne:
    settings['airborne'] = True
    
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


# grab calval data from json file
with open(cal_file,"r") as f:
    cal_json = json.loads(f.read())
f.close()

#default_aircraft_basepath = {
#    'CSET':'/scr/raf_data/CSET/24_Mar_17_BU/',
#    'SOCRATES':'/scr/raf_data/SOCRATES/' #SOCRATEStf01.nc       
#    } 
       
# Local Paths
aircraft_basepath = {
    'CSET':'/Users/mhayman/Documents/HSRL/aircraft_data/',
    'SOCRATES':'/Users/mhayman/Documents/HSRL/aircraft_data/' #SOCRATEStf01.nc       
    } 

if settings['airborne']:
    # list of aircraft variables to load
    var_aircraft = ['Time','GGALT','ROLL','PITCH','THDG','GGLAT','GGLON','TASX','ATX','PSXC']

    proj_list = []
    year_list = []
    for ai in range(len(cal_json['Flights'])):
        if not cal_json['Flights'][ai]['Project'] in proj_list:
            proj_list.extend([cal_json['Flights'][ai]['Project']])
            year_list.extend([lp.json_str_to_datetime(cal_json['Flights'][ai]['date'])])
            print('%d.) '%(len(proj_list)) + proj_list[-1] + ', ' + year_list[-1].strftime('%Y'))
    print('')
   
    # interactive prompt to determine desired flight
  
    usr_proj = np.int(input('Select Project: '))-1
    if usr_proj < 0 or usr_proj > len(proj_list)-1:
        print('Selection is not recognized')
    else:
        proj = proj_list[usr_proj]
    
    flight_list = []
    flight_date = []
    flight_label = []

                    
    for ai in range(len(cal_json['Flights'])):
        if cal_json['Flights'][ai]['Project'] == proj:
            flight_list.extend([proj+cal_json['Flights'][ai]['Flight Designation'] + str(cal_json['Flights'][ai]['Flight Number']).zfill(2)])
            flight_date.extend([lp.json_str_to_datetime(cal_json['Flights'][ai]['date'])])
            flight_label.extend([cal_json['Flights'][ai]['Flight Designation'].upper()+str(cal_json['Flights'][ai]['Flight Number']).zfill(2)])
            print('%d.) '%len(flight_list) + ' ' + flight_list[-1] + ', ' + flight_date[-1].strftime('%d-%b, %Y'))
    
    usr_flt = np.int(input('Select Flight: '))-1
    if usr_flt < 0 or usr_flt > len(flight_list)-1:
        print('Selection is not recognized')
    else:
        flt = flight_list[usr_flt]
            
    filePathAircraft = aircraft_basepath[proj] + flt + '.nc'
    
    
            
    #  load aircraft data    
    air_data = gv.load_aircraft_data(filePathAircraft,var_aircraft)
    
    # locate time range where aircraft is flying
    iflight = np.nonzero(air_data['TASX'] > settings['Airspeed_Threshold'])[0]
    it0 = iflight[0]  # index when aircraft starts moving
    it1 = iflight[-1]  # index when aircraft stops moving
    time_takeoff = flight_date[usr_flt]+datetime.timedelta(seconds=np.int(air_data['Time'][it0]))
    time_landing = flight_date[usr_flt]+datetime.timedelta(seconds=np.int(air_data['Time'][it1]))
    print('Flight time is: ')
    print('   '+time_takeoff.strftime('%H:%M %d-%b, %Y to'))
    print('   '+time_landing.strftime('%H:%M %d-%b, %Y'))
    print('')
    
    
    time_sel = input('Enter start hour (UTC) or press [Enter] select start of flight: ')
    if time_sel == '':
        time_start = flight_date[usr_flt]+datetime.timedelta(seconds=np.int(air_data['Time'][it0]))
        
    else:
        start_hr = np.float(time_sel)
        time_start = flight_date[usr_flt]+datetime.timedelta(seconds=np.int(start_hr*3600))
        
    
    time_sel = input('Enter duration or press [Enter] for end of flight: ')
    if time_sel == '':
        time_stop = flight_date[usr_flt]+datetime.timedelta(seconds=np.int(air_data['Time'][it1]))
    else:
        stop_hr = np.float(time_sel)
        time_stop = time_start+datetime.timedelta(seconds=np.int(stop_hr*3600))
    
    # check for out of bounds time limits
    if time_start > time_landing:
        time_start = time_landing - datetime.timedelta(minutes=1)
    if time_stop < time_takeoff:
        time_stop = time_takeoff + datetime.timedelta(minutes=1)
        
    cal_start = time_start
    cal_stop = time_stop

# grab raw data from netcdf files
[timeD,time_dt,time_sec],var_1d_data, profs = gv.load_raw_data(cal_start,cal_stop,var_2d_list,var_1d_list,basepath=basepath,verbose=True,as_prof=True)
#time_start = datetime.datetime(year_in,month_in,day_in)+datetime.timedelta(hours=start_hr)
#time_stop = time_start + datetime.timedelta(hours=stop_hr)


## grab raw data from netcdf files
#[timeD,time_dt,time_sec],var_1d_data, profs = gv.load_raw_data(time_start,time_stop,var_2d_list,var_1d_list,basepath=basepath,verbose=True,as_prof=True)

# find instances in raw data where I2 cell is removed
if 'RemoveLongI2Cell' in var_1d_data.keys():
    cal_indices = np.nonzero(var_1d_data['RemoveLongI2Cell'] < 50)[0]
else:
    cal_indices = []

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
#if settings['airborne']:
#    MaxAlt = np.nanmin(air_data['GGALT'])-settings['ground_range_buffer']
#else:
#    MaxAlt = settings['MaxAlt']
    
    
# set the master time to match all 2D profiles to
# (1d data will not be resampled)
master_time = np.arange(time_sec[0]-tres/2,time_sec[-1]+tres/2,tres)

time_1d,var_1d = gv.var_time_resample(master_time,time_sec,var_1d_data,average=True)
#var_1d = var_1d_data
#time_1d = time_sec

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
        
if settings['airborne']:
    air_data_t = gv.interp_aircraft_data(profs['combined_hi'].time,air_data)
    ground_mask = (air_data_t['GGALT']-settings['ground_range_buffer'])[:,np.newaxis] < profs['combined_hi'].range_array[np.newaxis,:]
    profs['combined_hi'].profile[np.nonzero(ground_mask)] = np.nan
    profs['combined_hi'].mask(ground_mask)
    
#    profs[var].slice_range(range_lim=[0,MaxAlt])
#    int_profs[var].bg_subtract(BGIndex)
#    int_profs[var].slice_range(range_lim=[0,MaxAlt])



if settings['airborne']:
    temp,pres = gv.get_TP_from_aircraft(air_data,profs['molecular'],telescope_direction=var_1d['TelescopeDirection'])
elif settings['load_reanalysis']:
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
if settings['EstimateExtinction'] and not settings['airborne']:
        
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
    Tatm = np.ones((1,profs['molecular'].profile.shape[1]))


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

geo_smooth = gv.savitzky_golay(geo_raw.profile.flatten(), settings['sg_win'], settings['sg_order'], deriv=0) 
#lp.plotprofiles([geo_raw],varplot=True,scale='log')

plt.figure()
plt.plot(geo_raw.profile.flatten())
plt.plot(geo_smooth)
plt.ylim([-2e-12,2e-12])


plt.show(block=False)

i_const = np.int(input('Make constant above index (e.g. 1000): '))
i_const_max = np.nonzero(geo_raw.SNR().flatten() < 0.2*geo_raw.SNR().flatten()[i_const])[0]
i_const_max = i_const_max[np.nonzero(i_const_max > i_const)[0][0]]

if settings['airborne']:
    geo_const = np.nanmean(geo_raw.profile[0,i_const-4:i_const])
    geo_fit = geo_smooth.copy()
    geo_fit[i_const:] = geo_const
else:
    geo_fit = gv.fit_high_range(geo_smooth,geo_raw.profile_variance,i_const,i_const_max)

plt.plot(geo_fit,'g--')
plt.show(block=False)

save_cal = input("Save Geo Calibrations[y/n]")
    
if save_cal == 'y' or save_cal == 'Y':
    write_data = np.zeros((geo_fit.size,2))
    write_data[:,0] = geo_raw.range_array
    write_data[:,1] = geo_fit
    
    header_str = 'profile data from ' + time_dt[0].strftime('%d-%b-%y %H:%M')+' -->'+time_dt[-1].strftime('%H:%M UTC\n') 
    header_str = header_str+ ' file created ' + datetime.datetime.now().strftime('%d-%b-%y %H:%M\n')
    header_str = header_str+ 'lidar ratio assumed %.1f below %.1f m'%(settings['LRassumed'],geo_raw.range_array[i_const])
    header_str = header_str+ 'Range  geo_correction'
    
    save_cal_file = 'geofile_default_'+time_dt[0].strftime('%Y%m%dT%H%M.geo')
    print('Saving to:\n  '+save_file_path+save_cal_file) 
    np.savetxt(save_file_path+save_cal_file,write_data,fmt='%.4e\t%.4e',header=header_str)
    
    save_cal_file_ez = 'geofile_GVHSRL'+time_dt[0].strftime('%Y%m%d')   
    print('Saving python vars to:\n  '+save_path_ez+save_cal_file_ez+'.npz') 
    
#    if hasattr(geo_raw.profile_variance,'mask'):
#        geo_var = geo_raw.profile_variance.flatten().data
#    else:
#        geo_var = geo_raw.profile_variance.flatten()
    
    np.savez(save_path_ez+save_cal_file_ez,start_time=cal_start,stop_time=cal_stop, \
            geo_mol=geo_fit,geo_mol_var=geo_raw.profile_variance.flatten().data,\
            geo_raw = geo_raw.profile.flatten().data,\
            sg_win = settings['sg_win'], sg_order = settings['sg_order'], \
            range_array=geo_raw.range_array, \
            TelescopeDirection = np.nanmean(var_1d_data['TelescopeDirection']), \
            airborne=settings['airborne'],Nprof=profs['moleculr'].NumProfList)
    




