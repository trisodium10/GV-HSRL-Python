# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 15:15:04 2017

@author: mhayman
"""

import sys
import os

#filename = inspect.getframeinfo(inspect.currentframe()).filename
filename = __file__


# add the path to GVHSRLlib manually
library_path = os.path.abspath(filename+'/../../libraries/')
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

#import ExternalDataFunctions as ex

import MLELidarProfileFunctions as mle
    


# input and raw_input are the same for python 3 and 2 respectively
# this makes it so input always accepts a string
try:
    input=raw_input
except NameError:
    pass

cal_file_path = os.path.abspath(filename+'/../../calibrations/cal_files/')+'/'
cal_file = cal_file_path + 'gv_calvals.json'

reanalysis_path = os.path.abspath(filename+'/../../../external_data/')+'/'

#save_file_path = '/Users/mhayman/Documents/Python/Lidar/'
#save_file_path = '/h/eol/mhayman/HSRL/hsrl_processing/hsrl_configuration/projDir/calfiles/'
default_settings = {
    'tres':0.5,  # resolution in time in seconds (0.5 sec) before altitude correction
    'tres_post':1*60, # resolution after altitude correction -  set to zero to not use
    'zres':7.5,  # altitude resolution in meters (7.5 m minimum)
    
    #mol_gain = 1.133915#1.0728915  # gain adjustment to molecular channel
    
    # index for where to treat the profile as background only
    'BGIndex': -100, # negative number provides an index from the end of the array
    'platform':'airborne', # 'ground' or 'airborne'.  If 'airborne' it needs an aircraft netcdf.
    'MaxAlt':10e3,
    'MinAlt':-30,
    
    'RemoveCals':True,  # don't include instances where the I2 cell is removed
                        # scan files are not included in the file search so they
                        # are removed anyway
    
    'Remove_Off_Data':True, # remove data where the lidar appears to be off
                            # currently only works if RemoveCals is True
    
    'diff_geo_correct':True,  # apply differential overlap correction
    
    'load_reanalysis':False, # load T and P reanalysis from NCEP/NCAR Model
    
    'plot_2D':True,   # pcolor plot the BSR and depolarization profiles
    'plot_date':True,  # plot results in date time format.  Otherwise plots as hour floats
    'save_plots':False, # save the plot data
    
    'save_data':False, # save data as netcdf
    
    'Estimate_Mol_Gain':True, # use statistics on BSR to estimate the molecular gain
    
    'hsrl_rb_adjust':True, # adjust for Rayleigh Brillouin Spectrum
    
    'Denoise_Mol':False, # run PTV denoising on molecular channel
    
    
    'Airspeed_Threshold':15, # threshold for determining start and end of the flight (in m/s)
    
    'loadQWP':'fixed',  # load 'fixed','rotating', or 'all' QWP data
    
    'as_altitude':True # process in altitude centered format or range centered format
    }
#sg_win = 11
#sg_order = 5

# check if any settings have been defined.  If not, define it as an empty dict.
try: settings
except NameError: settings = {}
    
# If a paremeter isn't supplied, use the default setting
for param in default_settings.keys():
    if not param in settings.keys():
        settings[param] = default_settings[param]    
    
tres = settings['tres']
tres_post = settings['tres_post']
zres = settings['zres']
BGIndex = settings['BGIndex']
MaxAlt = settings['MaxAlt']
MinAlt = settings['MinAlt']

#default_basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path
default_basepath = '/scr/rain1/rsfdata/projects/socrates/hsrl/raw/'

try:
    basepath
except NameError:
    basepath = default_basepath
#basepath = '/scr/eldora1/HSRL_data/'  # old path - still works with link from HSRL_data to /hsrl/raw/
#basepath = '/Users/mhayman/Documents/HSRL/GVHSRL_data/'  # local computer data path
        
default_aircraft_basepath = {
    'CSET':'/scr/raf_data/CSET/24_Mar_17_BU/',
    'SOCRATES':'/scr/raf_data/SOCRATES/' #SOCRATEStf01.nc       
    } 

try:
    aircraft_basepath
except NameError:
    aircraft_basepath = default_aircraft_basepath


# list of 1D variables to load
var_1d_list = ['total_energy','RemoveLongI2Cell'\
    ,'TelescopeDirection','TelescopeLocked','polarization','DATA_shot_count','builduptime']  # 'DATA_shot_count'

# list of 2D variables (profiles) to load
var_2d_list = ['molecular','combined_hi','combined_lo','cross']

# list of aircraft variables to load
var_aircraft = ['Time','GGALT','ROLL','PITCH','THDG','GGLAT','GGLON','TASX','ATX','PSXC']

# grab calval data from json file
with open(cal_file,"r") as f:
    cal_json = json.loads(f.read())
f.close()

proj_list = []
year_list = []
for ai in range(len(cal_json['Flights'])):
    if not cal_json['Flights'][ai]['Project'] in proj_list:
        proj_list.extend([cal_json['Flights'][ai]['Project']])
        year_list.extend([lp.json_str_to_datetime(cal_json['Flights'][ai]['date'])])
        print('%d.) '%(len(proj_list)) + proj_list[-1] + ', ' + year_list[-1].strftime('%Y'))
print('')
# check if the project/flight has been passed in 
# if not, ask the user for it    
try: 
    proj
except NameError:    
    # interactive prompt to determine desired flight
  
    usr_proj = np.int(input('Select Project: '))-1
    if usr_proj < 0 or usr_proj > len(proj_list)-1:
        print('Selection is not recognized')
    else:
        proj = proj_list[usr_proj]

flight_list = []
flight_date = []
flight_label = []
try:
    flt
    for ai in range(len(cal_json['Flights'])):
        if cal_json['Flights'][ai]['Project'] == proj:
            flight_list.extend([proj+cal_json['Flights'][ai]['Flight Designation'] + str(cal_json['Flights'][ai]['Flight Number']).zfill(2)])
            flight_date.extend([lp.json_str_to_datetime(cal_json['Flights'][ai]['date'])])
            flight_label.extend([cal_json['Flights'][ai]['Flight Designation'].upper()+str(cal_json['Flights'][ai]['Flight Number']).zfill(2)])
            print('%d.) '%len(flight_list) + ' ' + flight_list[-1] + ', ' + flight_date[-1].strftime('%d-%b, %Y'))
            if flight_list[-1] == flt:
                usr_flt = len(flight_list)-1
                print('^^^^^^^^^^')
                
except NameError:
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
print('Flight time is: ')
print('   '+(flight_date[usr_flt]+datetime.timedelta(seconds=np.int(air_data['Time'][it0]))).strftime('%H:%M %d-%b, %Y to'))
print('   '+(flight_date[usr_flt]+datetime.timedelta(seconds=np.int(air_data['Time'][it1]))).strftime('%H:%M %d-%b, %Y'))
print('')

try: 
    time_start = flight_date[usr_flt]+flight_time_start
    time_stop = flight_date[usr_flt]+flight_time_stop
except NameError:
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


print('Processing: ')
print('   '+time_start.strftime('%H:%M %d-%b, %Y to'))
print('   '+time_stop.strftime('%H:%M %d-%b, %Y'))
print('')

if settings['save_data']:
    try:
        save_data_file = save_data_path+flt+'_GVHSRL_'+time_start.strftime('%Y%m%dT%H%M')+'_'+time_stop.strftime('%Y%m%dT%H%M')+'.nc'
    except NameError:
        print('Save data is disabled')
        print('  No save path (save_data_path) is provided')
        settings['save_data'] = False

if settings['save_plots']:
    try:
        save_plots_path
        save_plots_base = flt+'_GVHSRL_'+time_start.strftime('%Y%m%dT%H%M')+'_'+time_stop.strftime('%Y%m%dT%H%M')
    except NameError:
        print('Save plots is disabled')
        print('  No save path (save_plots_path) is provided')
        settings['save_plots'] = False

# grab raw data from netcdf files
[timeD,time_dt,time_sec],var_1d_data, profs = gv.load_raw_data(time_start,time_stop,var_2d_list,var_1d_list,basepath=basepath,verbose=True,as_prof=True,loadQWP=settings['loadQWP'])


# find instances in raw data where I2 cell is removed
if 'RemoveLongI2Cell' in var_1d_data.keys():
    cal_indices = np.nonzero(var_1d_data['RemoveLongI2Cell'] < 50)[0]
else:
    cal_indices = []

if settings['Remove_Off_Data']:
    _,off_indices = profs['combined_hi'].trim_to_on(ret_index=True,delete=False)
    cal_indices = np.unique(np.concatenate((off_indices,cal_indices)))

# grab calibration data files
# down_gain is the molecular gain when the telescope points down
if settings['hsrl_rb_adjust']:
    mol_gain_up,diff_geo_file,mol_gain_down,diff_geo_file_down = lp.get_calval(time_start,cal_json,'Molecular Gain',cond=[['RB_Corrected','=','True']],returnlist=['value','diff_geo','down_gain','down_diff'])  
else:
    mol_gain_up,diff_geo_file,mol_gain_down,diff_geo_file_down = lp.get_calval(time_start,cal_json,"Molecular Gain",returnlist=['value','diff_geo','down_gain','down_diff'])
    
baseline_file = lp.get_calval(time_start,cal_json,"Baseline File")[0]
diff_pol_file = lp.get_calval(time_start,cal_json,"Polarization",returnlist=['diff_geo'])
i2_file = lp.get_calval(time_start,cal_json,"I2 Scan")

# load differential overlap correction
diff_data_up = np.load(cal_file_path+diff_geo_file)
if len(diff_geo_file_down) > 0:
    diff_data_down = np.load(cal_file_path+diff_geo_file_down)
else:
    diff_data = diff_data_up

baseline_data = np.load(cal_file_path+baseline_file)
if len(diff_pol_file):
    diff_pol_data = np.load(cal_file_path+diff_pol_file[0])

# load i2 scan from file
if len(i2_file):
    i2_data = np.load(cal_file_path+i2_file[0])
else:
    # if no i2 scan available, don't correct for Rayleigh Brillouin spectrum
    hsrl_rb_adjust = False

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

if settings['RemoveCals']:
    time_1d = np.delete(time_1d,cal_indices)
    var_1d = gv.delete_indices(var_1d,cal_indices)

if settings['as_altitude']:
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

# setup molecular gain vector based on telescope pointing direction
mol_gain = np.zeros(var_post['TelescopeDirection'].shape)
mol_gain[np.nonzero(var_post['TelescopeDirection']==1.0)] = mol_gain_up
mol_gain[np.nonzero(var_post['TelescopeDirection']==0.0)] = mol_gain_down
mol_gain = mol_gain[:,np.newaxis]

# setup variable diffierential overlap (up vs down pointing) if supplied
if len(diff_geo_file_down) > 0:
    diff_data = {}
    key_list = ['hi_diff_geo','lo_diff_geo']
    for var in key_list:
        diff_data[var] = np.ones((var_1d['TelescopeDirection'].size,diff_data_up[var].size))
        diff_data[var][np.nonzero(var_1d['TelescopeDirection']==1.0)[0],:] = diff_data_up[var]
        diff_data[var][np.nonzero(var_1d['TelescopeDirection']==0.0)[0],:] = diff_data_down[var]

"""
Main Profile Processing Loop
loop through each lidar profile in profs and perform basic processing 
operations
"""
    

int_profs = {}  # obtain time integrated profiles
for var in profs.keys():
    if settings['RemoveCals']:
        # remove instances where the I2 cell is removed
        profs[var].remove_time_indices(cal_indices)
#        profs[var].trim_to_on()  # remove points where lidar isn't transmitting
    if tres > 0.5:
        profs[var].time_resample(tedges=master_time,update=True,remainder=False)
    int_profs[var] = profs[var].copy()
    int_profs[var].time_integrate()
    
    if settings['as_altitude']:
        # maximum range required for this dataset
        range_trim = np.max([np.max(MaxAlt-air_data_t['GGALT']),np.max(air_data_t['GGALT']-MinAlt)])+4*zres
    else:
        # set maximum range to MaxAlt if range centered processing
        range_trim = MaxAlt
    
    if var == 'molecular' and settings['Denoise_Mol']:
        MolRaw = profs['molecular'].copy()    
    
    # background subtract the profile
    profs[var].bg_subtract(BGIndex)
    
    # profile specific processing routines
    if var == 'combined_hi' and settings['diff_geo_correct']:
        profs[var].diff_geo_overlap_correct(diff_data['hi_diff_geo'],geo_reference='molecular')
    elif var == 'combined_lo' and settings['diff_geo_correct']:
        profs[var].diff_geo_overlap_correct(diff_data['lo_diff_geo'],geo_reference='molecular')
        try:
            profs[var].gain_scale(1.0/diff_data['lo_norm'])
        except KeyError:
            print('No lo gain scale factor')
            print('   lo gain channel may need adjustment')
    elif var == 'cross' and settings['diff_geo_correct']:
        if len(diff_pol_file):
            profs[var].diff_geo_overlap_correct(diff_pol_data['cross_diff_geo'],geo_reference='combined_hi')
            profs[var].diff_geo_overlap_correct(diff_data['hi_diff_geo'],geo_reference='molecular')  
        else:
            profs[var].diff_geo_overlap_correct(diff_data['hi_diff_geo'],geo_reference='molecular')
    profs[var].slice_range(range_lim=[0,range_trim])
  
    if settings['as_altitude']:
        profs[var].range2alt(master_alt,air_data_t,telescope_direction=var_1d['TelescopeDirection'])
       
    
    if tres_post > 0 or tres <= 0.5:
        profs[var].time_resample(tedges=master_time_post,update=True,remainder=False)
        
    
    int_profs[var] = profs[var].copy()
    int_profs[var].time_integrate()



#if load_reanalysis:
#    pres,temp = ex.load_fixed_point_NCEP_TandP(profs['molecular'],lidar_location,reanalysis_path)


lp.plotprofiles(profs)
if settings['as_altitude']:
    temp,pres = gv.get_TP_from_aircraft(air_data,profs['molecular'])
else:
    temp,pres = gv.get_TP_from_aircraft(air_data,profs['molecular'],telescope_direction=var_post['TelescopeDirection'])
beta_m = lp.get_beta_m(temp,pres,profs['molecular'].wavelength)

if settings['Denoise_Mol']:
    MolDenoise,tune_list = mle.DenoiseMolecular(MolRaw,beta_m_sonde=beta_m, \
                            MaxAlt=MaxAlt,accel = False,tv_lim =[1.5, 2.8],N_tv_pts=59, \
                            geo_data=dict(geo_prof=np.array([2e14])),bg_index=-10,n=20)

if settings['hsrl_rb_adjust']:
    print('Obtaining Rayleigh-Brillouin Correction')
    dnu = 20e6  # resolution
    nu_max = 10e9 # max frequency relative to line center
    nu = np.arange(-nu_max,nu_max,dnu)
    Ti2 = np.interp(nu,i2_data['freq']*1e9,i2_data['mol_scan'])  # molecular transmission
    
    Tc2 = np.interp(nu,i2_data['freq']*1e9,i2_data['combined_scan'])  # combined transmission
    
    beta_mol_norm = lp.RB_Spectrum(temp.profile.flatten(),pres.profile.flatten()*9.86923e-6,profs['molecular'].wavelength,nu=nu,norm=True)
    eta_i2 = np.sum(Ti2[:,np.newaxis]*beta_mol_norm,axis=0)
    eta_i2 = eta_i2.reshape(temp.profile.shape)
    profs['molecular'].multiply_piecewise(1.0/eta_i2)

    profs['molecular'].gain_scale(mol_gain,gain_var = (mol_gain*0.05)**2)

        
    eta_c = np.sum(Tc2[:,np.newaxis]*beta_mol_norm,axis=0)
    eta_c = eta_c.reshape(temp.profile.shape)
    profs['combined_hi'].multiply_piecewise(1.0/eta_c)
    profs['cross'].multiply_piecewise(1.0/eta_c)
    
    if settings['Denoise_Mol']:
        MolDenoise.multiply_piecewise(1.0/eta_i2)
        MolDenoise.gain_scale(mol_gain)
else:
    # Rescale molecular channel to match combined channel gain
    profs['molecular'].gain_scale(mol_gain,gain_var = (mol_gain*0.05)**2)
    if settings['Denoise_Mol']:
        MolDenoise.gain_scale(mol_gain)


beta_a = lp.AerosolBackscatter(profs['molecular'],profs['combined_hi'],beta_m)



BSR = profs['combined_hi']/profs['molecular']
BSR.descript = 'Ratio of combined to molecular backscatter'
BSR.label = 'Backscatter Ratio'
BSR.profile_type = 'unitless'

#BSR_mask = (BSR.profile-1)/np.sqrt(BSR.profile_variance) < 5.0
#BSR_mask = BSR.profile < 2.0

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



if settings['Estimate_Mol_Gain']:
    # This segment estimates what the molecular gain should be 
    # based on a histogram minimum in BSR over the loaded data
    
    iUp = np.nonzero(var_post['TelescopeDirection']==1.0)
    
    BSRprof = BSR.profile[iUp,:].flatten()
    BSRalt = (np.ones(BSR.profile[iUp,:].shape)*BSR.range_array[np.newaxis,:]).flatten()
    
    BSRalt = np.delete(BSRalt,np.nonzero(np.isnan(BSRprof)))
    BSRprof = np.delete(BSRprof,np.nonzero(np.isnan(BSRprof)))

    bbsr = np.linspace(0,4,400)
    bsnr = np.linspace(1,100,250) # snr (70,100,250)

    balt = np.concatenate((BSR.range_array-BSR.mean_dR/2,BSR.range_array[-1:]+BSR.mean_dR/2))

    """
    perform analysis by altitude
    """
    
    hbsr = np.histogram2d(BSRprof,BSRalt,bins=[bbsr,balt])
    
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
    
    plt.figure()
    plt.pcolor(bbsr,balt,hbsr[0].T)
    plt.plot(hist_median,balt[1:],'r--')
    plt.plot(hist_med_sm,balt[1:],'g--')
    plt.xlabel('BSR')
    plt.ylabel('Altitude [m]')
    plt.title('Telescope Up')
    
    i_alt_lim = np.nonzero(np.logical_and(balt > 2000,balt < 6500))[0]
    
    mol_gain_adj = np.nanmin(hist_med_sm[i_alt_lim])
    
    print('Current Molecular (Telescope Up) Gain: %f'%mol_gain_up)
    print('Suggested Molecular (Telescope Up) Gain: %f'%(mol_gain_up*mol_gain_adj))
    
    iDown = np.nonzero(var_post['TelescopeDirection']==0.0)
    
    BSRprof = BSR.profile[iDown,:].flatten()
    BSRalt = (np.ones(BSR.profile[iDown,:].shape)*BSR.range_array[np.newaxis,:]).flatten()
    
    BSRalt = np.delete(BSRalt,np.nonzero(np.isnan(BSRprof)))
    BSRprof = np.delete(BSRprof,np.nonzero(np.isnan(BSRprof)))

    bbsr = np.linspace(0,4,400)
    bsnr = np.linspace(1,100,250) # snr (70,100,250)

    balt = np.concatenate((BSR.range_array-BSR.mean_dR/2,BSR.range_array[-1:]+BSR.mean_dR/2))

    """
    perform analysis by altitude
    """
    
    hbsr = np.histogram2d(BSRprof,BSRalt,bins=[bbsr,balt])
    
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
    
    plt.figure()
    plt.pcolor(bbsr,balt,hbsr[0].T)
    plt.plot(hist_median,balt[1:],'r--')
    plt.plot(hist_med_sm,balt[1:],'g--')
    plt.xlabel('BSR')
    plt.ylabel('Altitude [m]')
    plt.title('Telescope Down')
    
    i_alt_lim = np.nonzero(np.logical_and(balt > 2000,balt < 6500))[0]
    
    mol_gain_adj = np.nanmin(hist_med_sm[i_alt_lim])
    
    print('Current Molecular (Telescope Down) Gain: %f'%mol_gain_down)
    print('Suggested Molecular (Telescope Down) Gain: %f'%(mol_gain_down*mol_gain_adj))



# add a diagnostic for counts/backscatter coeff

# add a diagnostic for diff overlap between lo and hi channels as a function
# of count rate or backscatter coeff

#dPartMask = dPart.profile_variance > 1.0
#dPart.mask(dPartMask)
dPart.mask(dPart.profile_variance > 1.0)
dPart.mask(dPart.profile > 1.0)
dPart.mask(dPart.profile < -0.1)


proj_label = proj + ' ' + flight_label[usr_flt] + ', '

beta_a.mask(np.isnan(beta_a.profile))
BSR.mask(np.isnan(BSR.profile))
dPart.mask(np.isnan(dPart.profile))
dVol.mask(np.isnan(dVol.profile))

save_prof_list = [beta_a,dPart,dVol,BSR,beta_m]
save_var1d_post = {'TelescopeDirection':'1-Lidar Pointing Up, 0-Lidar Pointing Down',
                   'polarization':'System Quarter Waveplate orientation'}
save_air_post = {'THDG': 'aircraft heading',
                 'TASX': 'airspeed [m/s]',
                 'GGLAT': 'latitude',
                 'PITCH': 'pitch angle',
                 'GGALT': 'altitude [m]',
                 'PSXC': 'ambiant pressure in hPa',
                 'ROLL': 'aircraft roll angle',
                 'GGLON': 'longitude',
                 'ATX': 'ambiant temperature in C'}

if settings['save_data']:
    for ai in range(len(save_prof_list)):
        save_prof_list[ai].write2nc(save_data_file) #,name_override=True,tag=var_name)
    for var in save_var1d_post.keys():
        lp.write_var2nc(var_post[var],str(var),save_data_file,description=save_var1d_post[var])
    for var in save_air_post.keys():
        lp.write_var2nc(air_data_post[var],str(var),save_data_file,description=save_air_post[var])

if settings['plot_2D']:
    if settings['plot_date']:
        t1d_plt = x_time = mdates.date2num([datetime.datetime.fromordinal(BSR.StartDate.toordinal()) \
                    + datetime.timedelta(seconds=sec) for sec in time_1d])   
    else:
        t1d_plt = time_1d/3600.0
    
#    rfig = lp.pcolor_profiles([BSR,dPart],scale=['log','linear'],climits=[[1,1e2],[0,0.7]],ylimits=[MinAlt*1e-3,MaxAlt*1e-3],title_add=proj_label,plot_date=plot_date)
    rfig = lp.pcolor_profiles([beta_a,dPart],scale=['log','linear'],climits=[[1e-8,1e-3],[0,1.0]],ylimits=[MinAlt*1e-3,MaxAlt*1e-3],title_add=proj_label,plot_date=settings['plot_date'])
    for ai in range(len(rfig[1])):
        rfig[1][ai].plot(t1d_plt,air_data_t['GGALT']*1e-3,color='gray',linewidth=1.2)  # add aircraft altitude
        
    if settings['save_plots']:
        plt.savefig(save_plots_path+'Aerosol_Backscatter_'+save_plots_base,dpi=300)
        
    rfig = lp.pcolor_profiles([profs['combined_hi'],dVol],scale=['log','linear'],climits=[[1e-1,1e4],[0,1.0]],ylimits=[MinAlt*1e-3,MaxAlt*1e-3],title_add=proj_label,plot_date=settings['plot_date'])
    for ai in range(len(rfig[1])):
        rfig[1][ai].plot(t1d_plt,air_data_t['GGALT']*1e-3,color='gray',linewidth=1.2)  # add aircraft altitude
    if settings['save_plots']:
        plt.savefig(save_plots_path+'AttenuatedBackscatter_'+save_plots_base,dpi=300)
    #lp.plotprofiles(profs)
    #dPart.mask(dPartMask)
    #lp.pcolor_profiles([BSR,dVol],scale=['log','linear'],climits=[[1,5e2],[0,1.0]])
    #lp.pcolor_profiles([dVol],scale=['linear'],climits=[[0,1]])
plt.show()