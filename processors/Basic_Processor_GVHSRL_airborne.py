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

import matplotlib.dates as mdates

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

tres = 0.5  # resolution in time in seconds (0.5 sec) before altitude correction
tres_post = 1*60 # resolution after altitude correction -  set to zero to not use
zres = 7.5  # altitude resolution in meters (7.5 m minimum)

#mol_gain = 1.133915#1.0728915  # gain adjustment to molecular channel

# index for where to treat the profile as background only
BGIndex = -100; # negative number provides an index from the end of the array
platform = 'ground' # 'ground' or 'airborne'.  If 'airborne' it needs an aircraft netcdf.
MaxAlt = 10e3
MinAlt = -30

RemoveCals = True   # don't include instances where the I2 cell is removed
                    # scan files are not included in the file search so they
                    # are removed anyway

diff_geo_correct = False  # apply differential overlap correction

load_reanalysis = False # load T and P reanalysis from NCEP/NCAR Model

plot_2D = True   # pcolor plot the BSR and depolarization profiles
plot_date = True  # plot results in date time format.  Otherwise plots as hour floats

Estimate_Mol_Gain = True # use statistics on BSR to estimate the molecular gain


Airspeed_Threshold = 15 # threshold for determining start and end of the flight (in m/s)

#sg_win = 11
#sg_order = 5

#basepath = '/scr/eldora1/HSRL_data/'  # old path - still works with link from HSRL_data to /hsrl/raw/
basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path
#basepath = '/Users/mhayman/Documents/HSRL/GVHSRL_data/'  # local computer data path
        
aircraft_basepath = '/scr/raf_data/CSET/24_Mar_17_BU/'
        



# list of 1D variables to load
var_1d_list = ['total_energy','RemoveLongI2Cell'\
    ,'TelescopeDirection','TelescopeLocked','polarization','DATA_shot_count']  # 'DATA_shot_count'

# list of 2D variables (profiles) to load
var_2d_list = ['molecular','combined_hi','combined_lo','cross']

# list of aircraft variables to load
var_aircraft = ['Time','GGALT','ROLL','PITCH','THDG','GGLAT','GGLON','TASX','ATX','PSXC']

# grab calval data from json file
with open(cal_file,"r") as f:
    cal_json = json.loads(f.read())
f.close()

# interactive prompt to determine desired flight

proj_list = []
year_list = []
for ai in range(len(cal_json['Flights'])):
    if not cal_json['Flights'][ai]['Project'] in proj_list:
        proj_list.extend([cal_json['Flights'][ai]['Project']])
        year_list.extend([lp.json_str_to_datetime(cal_json['Flights'][ai]['date'])])
        print('%d.) '%(len(proj_list)) + proj_list[-1] + ', ' + year_list[-1].strftime('%Y'))

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
    filePathAircraft = aircraft_basepath + flight_list[usr_flt] + '.nc'
    
#  load aircraft data
    
air_data = gv.load_aircraft_data(filePathAircraft,var_aircraft)

# locate time range where aircraft is flying
iflight = np.nonzero(air_data['TASX'] > Airspeed_Threshold)[0]
it0 = iflight[0]  # index when aircraft starts moving
it1 = iflight[-1]  # index when aircraft stops moving
print('Flight time is: ')
print('   '+(flight_date[usr_flt]+datetime.timedelta(seconds=np.int(air_data['Time'][it0]))).strftime('%H:%M %d-%b, %Y to'))
print('   '+(flight_date[usr_flt]+datetime.timedelta(seconds=np.int(air_data['Time'][it1]))).strftime('%H:%M %d-%b, %Y'))
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




# grab raw data from netcdf files
[timeD,time_dt,time_sec],var_1d_data, profs = gv.load_raw_data(time_start,time_stop,var_2d_list,var_1d_list,basepath=basepath,verbose=True,as_prof=True)


# find instances in raw data where I2 cell is removed
if 'RemoveLongI2Cell' in var_1d_data.keys():
    cal_indices = np.nonzero(var_1d_data['RemoveLongI2Cell'] < 50)[0]
else:
    cal_indices = []



mol_gain,diff_geo_file = lp.get_calval(time_start,cal_json,"Molecular Gain",returnlist=['value','diff_geo'])
baseline_file = lp.get_calval(time_start,cal_json,"Baseline File")[0]
# load differential overlap correction
diff_data = np.load(cal_file_path+diff_geo_file)
baseline_data = np.load(cal_file_path+baseline_file)
#diff_data = np.load(cal_file_path+'diff_geo_GVHSRL20171025_tmp.npz')
#baseline_data = np.load(cal_file_path+'diff_geo_GVHSRL20171025_tmp.npz')

lidar_location = lp.get_calval(time_start,cal_json,"Location",returnlist=['latitude','longitude'])

# set the master time to match all 2D profiles to
# (1d data will not be resampled)
#master_time = np.arange(time_sec[0]-tres/2,time_sec[-1]+tres/2,tres) #
sec_start = np.max([time_sec[0],(time_start-flight_date[usr_flt]).total_seconds()])
sec_stop = np.min([time_sec[-1],(time_stop-flight_date[usr_flt]).total_seconds()])
if tres > 0.5:
    master_time = np.arange(sec_start-tres/2,sec_stop+tres/2,tres)
    time_1d,var_1d = gv.var_time_resample(master_time,time_sec,var_1d_data,average=True)
else:
    time_1d = time_sec.copy()
    var_1d = var_1d_data.copy()

master_alt = np.arange(MinAlt,MaxAlt+zres,zres)


air_data_t = gv.interp_aircraft_data(time_1d,air_data)

# time resolution after range to altitude conversion
if tres_post > 0:
    master_time_post = np.arange(sec_start-tres_post/2,sec_stop+tres_post/2,tres_post)
    time_post,var_post = gv.var_time_resample(master_time_post,time_sec,var_1d_data,average=True)
    air_data_post = gv.interp_aircraft_data(time_post,air_data)
elif tres > 0.5:
    master_time_post = master_time
    time_post = time_1d
    var_post = var_1d
    air_data_post = air_data_t
else:
    master_time_post = np.arange(sec_start-tres/2,sec_stop+tres/2,tres)
    

int_profs = {}  # obtain time integrated profiles
for var in profs.keys():
    if RemoveCals:
        # remove instances where the I2 cell is removed
        profs[var].remove_time_indices(cal_indices)
    if tres > 0.5:
        profs[var].time_resample(tedges=master_time,update=True,remainder=False)
    int_profs[var] = profs[var].copy()
    int_profs[var].time_integrate()
    
    # maximum range required for this dataset
    range_trim = np.max([np.max(MaxAlt-air_data_t['GGALT']),np.max(air_data_t['GGALT']-MinAlt)])+4*zres
    
    profs[var].bg_subtract(BGIndex)
    if var == 'combined_hi' and diff_geo_correct:
        profs[var].diff_geo_overlap_correct(diff_data['hi_diff_geo'])
    elif var == 'combined_lo' and diff_geo_correct:
        profs[var].diff_geo_overlap_correct(diff_data['hi_diff_geo'])
        profs[var].gain_scale(1.0/diff_data['lo_norm'])
    profs[var].slice_range(range_lim=[0,range_trim])
    profs[var].range2alt(master_alt,air_data_t,telescope_direction=var_1d['TelescopeDirection'])
    if tres_post > 0 or tres <= 0.5:
        profs[var].time_resample(tedges=master_time_post,update=True,remainder=False)
    
    int_profs[var] = profs[var].copy()
    int_profs[var].time_integrate()
#    int_profs[var].bg_subtract(BGIndex)
#    int_profs[var].slice_range(range_lim=[0,MaxAlt])


if load_reanalysis:
    pres,temp = ex.load_fixed_point_NCEP_TandP(profs['molecular'],lidar_location,reanalysis_path)


lp.plotprofiles(profs)

Temp,Pres = gv.get_TP_from_aircraft(air_data,profs['molecular'])
beta_m = lp.get_beta_m(Temp,Pres,profs['molecular'].wavelength)


profs['molecular'].gain_scale(mol_gain,gain_var = (mol_gain*0.1)**2)
#profs['combined_hi'].diff_geo_overlap_correct(diff_data['hi_diff_geo'][:profs['combined_hi'].range_array.size])
#profs['combined_lo'].diff_geo_overlap_correct(diff_data['lo_diff_geo'][:profs['combined_lo'].range_array.size])
#profs['combined_lo'].gain_scale(1.0/diff_data['lo_norm'])


beta_a = lp.AerosolBackscatter(profs['molecular'],profs['combined_hi'],beta_m)

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
    bsnr = np.linspace(10,380,100)
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

dPartMask = dPart.profile_variance > 1.0
dPart.mask(dPartMask)

proj_label = proj + ' ' + flight_label[usr_flt] + ', '

BSR.mask(np.isnan(BSR.profile))
dPart.mask(np.isnan(dPart.profile))
dVol.mask(np.isnan(dVol.profile))

if plot_2D:
    if plot_date:
        t1d_plt = x_time = mdates.date2num([datetime.datetime.fromordinal(BSR.StartDate.toordinal()) \
                    + datetime.timedelta(seconds=sec) for sec in time_1d])   
    else:
        t1d_plt = time_1d/3600.0
    
#    rfig = lp.pcolor_profiles([BSR,dPart],scale=['log','linear'],climits=[[1,1e2],[0,0.7]],ylimits=[MinAlt*1e-3,MaxAlt*1e-3],title_add=proj_label,plot_date=plot_date)
    rfig = lp.pcolor_profiles([beta_a,dPart],scale=['log','linear'],climits=[[1e-8,1e-3],[0,0.7]],ylimits=[MinAlt*1e-3,MaxAlt*1e-3],title_add=proj_label,plot_date=plot_date)
    for ai in range(len(rfig[1])):
        rfig[1][ai].plot(t1d_plt,air_data_t['GGALT']*1e-3,color='gray',linewidth=1.2)  # add aircraft altitude
        
    rfig = lp.pcolor_profiles([profs['combined_hi'],dVol],scale=['log','linear'],climits=[[1e-1,1e4],[0,1.0]],ylimits=[MinAlt*1e-3,MaxAlt*1e-3],title_add=proj_label,plot_date=plot_date)
    for ai in range(len(rfig[1])):
        rfig[1][ai].plot(t1d_plt,air_data_t['GGALT']*1e-3,color='gray',linewidth=1.2)  # add aircraft altitude
    #lp.plotprofiles(profs)
    #dPart.mask(dPartMask)
    #lp.pcolor_profiles([BSR,dVol],scale=['log','linear'],climits=[[1,5e2],[0,1.0]])
    #lp.pcolor_profiles([dVol],scale=['linear'],climits=[[0,1]])
