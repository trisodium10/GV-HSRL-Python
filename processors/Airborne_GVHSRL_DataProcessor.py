# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 13:59:17 2018

@author: mhayman
"""

#import sys
import os

#filename = inspect.getframeinfo(inspect.currentframe()).filename
#try:
#    filename = fullpath
#except NameError:
#    filename = __file__


## add the path to GVHSRLlib manually
#library_path = os.path.abspath(filename+'/../../libraries/')
#print(library_path)
#if library_path not in sys.path:
#    sys.path.append(library_path)


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

def ProcessAirborneDataChunk(time_start,time_stop,
                             settings={},paths={},process_vars={},
                             date_reference=0):

    filename = __file__
#    print(filename)
    cal_file_path = os.path.abspath(filename+'/../../calibrations/cal_files/')+'/'
    cal_file = cal_file_path + 'gv_calvals.json'
    
#    reanalysis_path = os.path.abspath(filename+'/../../../external_data/')+'/'
    
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
        
        'get_extinction':False,  # process data for extinction   
        'ext_sg_width':21,  # savitsky-gouley window width (must be odd)
        'ext_sg_order':3,   # savitsky-gouley polynomial order
        'ext_tres':15,      # extinction convolution kernel time width in seconds
        'ext_zres':60,      # extinction convolution kernel range width in meters
        
        'baseline_subtract':False, # use baseline subtraction
        'deadtime_correct':False,  # correct for APD deadtime
        
        'time_ref2takeoff':False,    # flight_time_start and 
                                     # time_stop are referened to takeoff time
        
        'RemoveCals':True,  # don't include instances where the I2 cell is removed
                            # scan files are not included in the file search so they
                            # are removed anyway
        
        'Remove_Off_Data':True, # remove data where the lidar appears to be off
                                # currently only works if RemoveCals is True
        
        'diff_geo_correct':True,  # apply differential overlap correction
        
        'load_reanalysis':False, # load T and P reanalysis from NCEP/NCAR Model
        
        'plot_2D':True,   # pcolor plot the BSR and depolarization profiles
        'show_plots':True, # show plots in matplotlib window
        'plot_date':True,  # plot results in date time format.  Otherwise plots as hour floats
        'save_plots':False, # save the plot data
        
        'save_data':False, # save data as netcdf
        
        'time_axis_scale':5.0,  # scale for horizontal axis on pcolor plots
        'alt_axis_scale':1.0,   # scale for vertical axis on pcolor plots
        'count_mask_threshold':2.0,  # count mask threshold (combined_hi).  If set to zero, no mask applied  
        'd_part_res_lim':0.25,  # resolution limit to decide where to mask particle depolarization data
        
        'Estimate_Mol_Gain':True, # use statistics on BSR to estimate the molecular gain
        
        'hsrl_rb_adjust':True, # adjust for Rayleigh Brillouin Spectrum
        
        'Denoise_Mol':False, # run PTV denoising on molecular channel
        
        
        'Airspeed_Threshold':15, # threshold for determining start and end of the flight (in m/s)
        
        'loadQWP':'fixed',  # load 'fixed','rotating', or 'all' QWP data
        
        'as_altitude':False, # process in altitude centered format or range centered format
        
        'SNRlimit':40.0,  # minimum integrated SNR to treat the lidar as transmitting
                         # used to filter instances where the shutter is closed
                         # toggle this with 'Remove_Off_Data'
        
        'use_aircraft_tref':True  # set the time reference based on aircraft data
        }
    
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
    
    if not hasattr(date_reference,'year'):
        date_reference = datetime.datetime(year=time_start.year,month=time_start.month,day=time_start.day)
        
    
    #default_basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path
    #default_basepath = '/scr/rain1/rsfdata/projects/socrates/hsrl/raw/'
    default_basepath = '/Users/mhayman/Documents/HSRL/GVHSRL_data/'  # local path
    
    try:
        basepath = paths['basepath']
    except KeyError:
        basepath = default_basepath
        paths['basepath'] = basepath
    #basepath = '/scr/eldora1/HSRL_data/'  # old path - still works with link from HSRL_data to /hsrl/raw/
    #basepath = '/Users/mhayman/Documents/HSRL/GVHSRL_data/'  # local computer data path
     
           
    #default_aircraft_basepath = {
    #    'CSET':'/scr/raf_data/CSET/24_Mar_17_BU/',
    #    'SOCRATES':'/scr/raf_data/SOCRATES/' #SOCRATEStf01.nc       
    #    } 
           
#    # Local Paths
#    default_aircraft_basepath = {
#        'CSET':'/Users/mhayman/Documents/HSRL/aircraft_data/',
#        'SOCRATES':'/Users/mhayman/Documents/HSRL/aircraft_data/' #SOCRATEStf01.nc       
#        } 
    
#    if 'aircraft_basepath' in paths.keys():
#        aircraft_basepath = paths['aircraft_basepath']
#    else:
#        aircraft_basepath = default_aircraft_basepath
    
    """
    Load variable lists
    """
    
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
    
       
    filePathAircraft = paths['filePathAircraft']
    #  load aircraft data    
    air_data,aircraft_t_ref = gv.load_aircraft_data(filePathAircraft,var_aircraft)

    print('Processing: ')
    print('   '+time_start.strftime('%H:%M %d-%b, %Y to'))
    print('   '+time_stop.strftime('%H:%M %d-%b, %Y'))
    print('')
    
    if settings['use_aircraft_tref']:
        date_reference = aircraft_t_ref
    
#    plt.figure()
#    plt.plot(air_data['Time']/3600.0,air_data['TASX'])
#    plt.show()    
    
    flt = process_vars['flt']
    
    if settings['save_data']:
        try:
            save_data_path = paths['save_data_path']
            save_data_file = save_data_path+flt+'_GVHSRL_'+time_start.strftime('%Y%m%dT%H%M')+'_'+time_stop.strftime('%Y%m%dT%H%M')+'.nc'
        except KeyError:
            print('Save data is disabled')
            print('  No save path (save_data_path) is provided')
            settings['save_data'] = False
    
    if settings['save_plots']:
        try:
            save_plots_path = paths['save_plots_path']
            save_plots_base = flt+'_GVHSRL_'+time_start.strftime('%Y%m%dT%H%M')+'_'+time_stop.strftime('%Y%m%dT%H%M')
        except KeyError:
            print('Save plots is disabled')
            print('  No save path (save_plots_path) is provided')
            settings['save_plots'] = False
    
    # grab raw data from netcdf files
    time_list,var_1d_data, profs = gv.load_raw_data(time_start,time_stop,var_2d_list,var_1d_list,basepath=basepath,verbose=True,as_prof=True,loadQWP=settings['loadQWP'],date_reference=date_reference)
    
    
    
    run_processing = len(profs) > 0    
    
    while run_processing:
        #execute processing if data was found
    
        timeD = time_list[0]
        time_dt = time_list[1]
        time_sec = time_list[2]
    
        # find instances in raw data where I2 cell is removed
        if 'RemoveLongI2Cell' in var_1d_data.keys():
            cal_indices = np.nonzero(var_1d_data['RemoveLongI2Cell'] < 50)[0]
        else:
            cal_indices = []
        
        # find instances where the lidar is not transmitting
        if settings['Remove_Off_Data']:
            _,off_indices = profs['combined_hi'].trim_to_on(ret_index=True,delete=False,SNRlim=settings['SNRlimit'])
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
        dead_time_list = lp.get_calval(time_start,cal_json,"Dead_Time",returnlist=['combined_hi','cross','combined_lo','molecular'])
        dead_time = dict(zip(['combined_hi','cross','combined_lo','molecular'],dead_time_list))
        
        if settings['get_extinction']:
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
        
        baseline_data = np.load(cal_file_path+baseline_file)
        
        if len(diff_pol_file):
            diff_pol_data = np.load(cal_file_path+diff_pol_file[0])
        
        # load i2 scan from file
        if len(i2_file):
            i2_data = np.load(cal_file_path+i2_file[0])
        else:
            # if no i2 scan available, don't correct for Rayleigh Brillouin spectrum
            settings['hsrl_rb_adjust'] = False
        
        lidar_location = lp.get_calval(time_start,cal_json,"Location",returnlist=['latitude','longitude'])
        
        
        flight_date = process_vars['flight_date']
        usr_flt = process_vars['usr_flt']
        
        print('flight_date: '+flight_date[usr_flt].strftime('%Y-%b-%d %H:%M'))
        print('date_reference: '+date_reference.strftime('%Y-%b-%d %H:%M'))
        
        # set the master time to match all 2D profiles to
        # (1d data will not be resampled)
        #master_time = np.arange(time_sec[0]-tres/2,time_sec[-1]+tres/2,tres) #
#        sec_start = np.max([time_sec[0],(time_start-flight_date[usr_flt]).total_seconds()])
#        sec_stop = np.min([time_sec[-1],(time_stop-flight_date[usr_flt]).total_seconds()])
#        sec_start = np.max([time_sec[0],(time_start-date_reference).total_seconds()])
#        sec_stop = np.min([time_sec[-1],(time_stop-date_reference).total_seconds()])
        sec_start = (time_start-date_reference).total_seconds()
        sec_stop = (time_stop-date_reference).total_seconds()
        print('found data for')
        print('   %f h-UTC to'%(sec_start/3600.0))
        print('   %f h-UTC to'%(sec_stop/3600.0))
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
            print('post time res: %f seconds'%tres_post)
            master_time_post = np.arange(sec_start-tres_post/2,sec_stop+tres_post/2,tres_post)
        #    time_post,var_post = gv.var_time_resample(master_time_post,time_sec,var_1d_data,average=True)
        #    air_data_post = gv.interp_aircraft_data(time_post,air_data)
        elif tres > 0.5:
            print('using tres')
            master_time_post = master_time
            time_post = time_1d
            var_post = var_1d
            air_data_post = air_data_t
        else:
            print('no change to time resolution')
            master_time_post = np.arange(sec_start-tres/2,sec_stop+tres/2,tres)
        print('master_time_post limits')
        print('   %f - %f h UTC'%(master_time_post[0]/3600.0,master_time_post[-1]/3600.0))
        print('   %f second resolution'%np.mean(np.diff(master_time_post)))
        print('   %d data points'%master_time_post.size)
        
        ## setup molecular gain vector based on telescope pointing direction
        #mol_gain = np.zeros(var_post['TelescopeDirection'].shape)
        #mol_gain[np.nonzero(var_post['TelescopeDirection']==1.0)] = mol_gain_up
        #mol_gain[np.nonzero(var_post['TelescopeDirection']==0.0)] = mol_gain_down
        #mol_gain = mol_gain[:,np.newaxis]
        
        # setup variable diffierential overlap (up vs down pointing) if supplied
        if len(diff_geo_file_down) > 0:
            diff_data = {}
            key_list = ['hi_diff_geo','lo_diff_geo']
            for var in key_list:
                diff_data[var] = np.ones((var_1d['TelescopeDirection'].size,diff_data_up[var].size))
                diff_data[var][np.nonzero(var_1d['TelescopeDirection']==1.0)[0],:] = diff_data_up[var]
                diff_data[var][np.nonzero(var_1d['TelescopeDirection']==0.0)[0],:] = diff_data_down[var]
        
        # setup variable geo overlap (up vs down pointing) if supplied
        if settings['get_extinction']:
            if len(geo_file_down) > 0:
                geo_data = {}
                key_list = ['geo_mol','geo_mol_var','Nprof']
                for var in key_list:
                    if var in geo_up.keys():
                        geo_data[var] = np.ones((var_1d['TelescopeDirection'].size,geo_up[var].size))
                        geo_data[var][np.nonzero(var_1d['TelescopeDirection']==1.0)[0],:] = geo_up[var]
                        if var in geo_down.keys():
                            geo_data[var][np.nonzero(var_1d['TelescopeDirection']==0.0)[0],:] = geo_down[var]
                        else:
                            geo_data[var][np.nonzero(var_1d['TelescopeDirection']==0.0)[0],:] = geo_up[var]
                    else:
                        geo_data[var] = np.ones((var_1d['TelescopeDirection'].size,1))
        
        """
        Main Profile Processing Loop
        loop through each lidar profile in profs and perform basic processing 
        operations
        """
        
        if settings['as_altitude']:
            # maximum range required for this dataset
            range_trim = np.max([np.max(MaxAlt-air_data_t['GGALT']),np.max(air_data_t['GGALT']-MinAlt)])+4*zres
        else:
            # set maximum range to MaxAlt if range centered processing
            range_trim = MaxAlt
        
        int_profs = {}  # obtain time integrated profiles
        for var in profs.keys():
            if settings['RemoveCals']:
                # remove instances where the I2 cell is removed
                profs[var].remove_time_indices(cal_indices)
        #        profs[var].trim_to_on()  # remove points where lidar isn't transmitting
            if tres > 0.5:
                profs[var].time_resample(tedges=master_time,update=True,remainder=False)
#            int_profs[var] = profs[var].copy()
#            int_profs[var].time_integrate()
            
            if settings['deadtime_correct']:
#                if var == 'combined_hi':
#                    p_before = profs[var].copy()
                if hasattr(profs[var],'NumProfsList'):
                    profs[var].nonlinear_correct(dead_time[var],laser_shot_count=2000*profs[var].NumProfsList[:,np.newaxis],std_deadtime=5e-9)
                else:
                    # number of laser shots is based on an assumption that there is one 0.5 second profile per time bin
                    profs[var].nonlinear_correct(dead_time[var],laser_shot_count=2000,std_deadtime=5e-9)

#                if var == 'combined_hi':
#                    p_after = profs[var].copy()
#                    lp.plotprofiles([p_before,p_after],varplot=True,time=18.3*3600)
#                    lp.plotprofiles([p_after],varplot=True,time=18.1*3600)
            
            if var == 'molecular' and settings['Denoise_Mol']:
                MolRaw = profs['molecular'].copy()    
    
            if settings['baseline_subtract']:        
                # baseline subtract  profiles
                profs[var].baseline_subtract(baseline_data['save_data'][var]['fit'], \
                    baseline_var = baseline_data['save_data'][var]['variance'], \
                    tx_norm=var_1d['total_energy'][:,np.newaxis]/baseline_data['avg_energy'])
                    
            
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
            elif var == 'molecular' and settings['get_extinction']:
                # if retrieving extinction, use a range centered profile to obtain it
                mol_ext = profs['molecular'].copy()
                mol_ext.multiply_piecewise(geo_data['geo_mol'])
                mol_ext.slice_range(range_lim=[0,range_trim])
                if tres_post > 0 or tres <= 0.5:
                    mol_ext.time_resample(tedges=master_time_post,update=True,remainder=False)
                
                
            profs[var].slice_range(range_lim=[0,range_trim])
          
            if settings['as_altitude']:
                profs[var].range2alt(master_alt,air_data_t,telescope_direction=var_1d['TelescopeDirection'])
               
            
            if tres_post > 0 or tres <= 0.5:
#                print(var+' resample:')
#                print('  [start,end] = [%f,%f] h UTC'%(profs[var].time[0]/3600.0,profs[var].time[-1]/3600.0))
#                print('  resolution = %f sec'%profs[var].mean_dt)
#                print('  %d data points'%profs[var].profile.shape[0])
                profs[var].time_resample(tedges=master_time_post,update=True,remainder=False)
            if profs[var].profile.size == 0:
                run_processing = False
#                print('  -- Now --')
#                print('  [start,end] = [%f,%f]'%(profs[var].time[0]/3600.0,profs[var].time[-1]/3600.0))
#                print('  resolution = %f sec'%profs[var].mean_dt)
#                print('  %d data points'%profs[var].profile.shape[0])
#                print(profs[var].ProcessingStatus)
                
            else:
                int_profs[var] = profs[var].copy()
                int_profs[var].time_integrate()
        if not run_processing:
            break
        # reformulate the master time based on the time that appears in the profiles
        # after processing
        if tres_post > 0:
#            print(master_time_post.shape)
            
            master_time_post = np.concatenate((np.array([profs['molecular'].time[0]-tres_post*0.5]), \
                0.5*np.diff(profs['molecular'].time)+profs['molecular'].time[:-1], \
                np.array([profs['molecular'].time[-1]+tres_post*0.5])))
#            print(master_time_post.shape)
            time_post,var_post = gv.var_time_resample(master_time_post,time_sec,var_1d_data,average=True)
            air_data_post = gv.interp_aircraft_data(time_post,air_data)
        else:
            time_post = time_1d
            var_post = var_1d
            air_data_post = air_data_t
        
        # setup molecular gain vector based on telescope pointing direction
        mol_gain = np.zeros(var_post['TelescopeDirection'].shape)
        mol_gain[np.nonzero(var_post['TelescopeDirection']==1.0)] = mol_gain_up
        mol_gain[np.nonzero(var_post['TelescopeDirection']==0.0)] = mol_gain_down
        mol_gain = mol_gain[:,np.newaxis]
        
        #if load_reanalysis:
        #    pres,temp = ex.load_fixed_point_NCEP_TandP(profs['molecular'],lidar_location,reanalysis_path)
        
        
        lp.plotprofiles(profs)
        if settings['as_altitude']:
            temp,pres = gv.get_TP_from_aircraft(air_data,profs['molecular'])
        else:
            temp,pres = gv.get_TP_from_aircraft(air_data,profs['molecular'],telescope_direction=var_post['TelescopeDirection'])
        beta_m = lp.get_beta_m(temp,pres,profs['molecular'].wavelength)
        
        
        if (settings['get_extinction']  or settings['Denoise_Mol']) and settings['as_altitude']:
#            print(var_post['TelescopeDirection'].size)
#            print(mol_ext.profile.shape)
#            print(profs['molecular'].time.size)
#            print(time_post.size)
    #        temp_ext,pres_ext = gv.get_TP_from_aircraft(air_data,mol_ext,telescope_direction=var_1d['TelescopeDirection'])
            temp_ext,pres_ext = gv.get_TP_from_aircraft(air_data,mol_ext,telescope_direction=var_post['TelescopeDirection'])
            beta_m_ext = lp.get_beta_m(temp_ext,pres_ext,profs['molecular'].wavelength)
        else:
            beta_m_ext = beta_m.copy()
        
        
        if settings['Denoise_Mol']:
            
            if len(geo_file_down) > 0:
                geo_denoise = {}
                key_list = ['geo_mol','geo_mol_var','Nprof']
                for var in key_list:
                    if var in geo_up.keys():
                        geo_denoise[var] = np.ones((var_1d['TelescopeDirection'].size,geo_up[var].size))
                        geo_denoise[var][np.nonzero(var_1d['TelescopeDirection']==1.0)[0],:] = geo_up[var]
                        if var in geo_down.keys():
                            geo_denoise[var][np.nonzero(var_1d['TelescopeDirection']==0.0)[0],:] = geo_down[var]
                        else:
                            geo_denoise[var][np.nonzero(var_1d['TelescopeDirection']==0.0)[0],:] = geo_up[var]
                    else:
                        geo_data[var] = np.ones((var_1d['TelescopeDirection'].size,1))            
                geo_denoise['range_array'] = geo_up['range_array'].copy()
                
            print('Denoising Molecular Channel')
            MolDenoise,tune_list = gv.DenoiseMolecular(MolRaw,beta_m_sonde=beta_m_ext.copy(), \
                                    MaxAlt=range_trim,accel = True,tv_lim =[1.5, 2.8],N_tv_pts=59, \
                                    bg_index=-10,n=1,geo_data=geo_denoise,geo_key='geo_mol',verbose=False) # dict(geo_prof=np.array([2e14])), geo_data=geo_data,geo_key='geo_mol'
#            # testing and debugging
#            MolRaw.bg_subtract(-10)
#            lp.plotprofiles([MolRaw,MolDenoise],time=22.1*3600)
##            lp.plotprofiles([MolRaw,MolDenoise],time=22.2*3600)
#            plt.show()
            
            MolDenoise.slice_range(range_lim=[0,range_trim])

            if settings['as_altitude']:
                MolDenoise.range2alt(master_alt,air_data_t,telescope_direction=var_1d['TelescopeDirection'])
            if tres_post > 0 or tres <= 0.5:
                MolDenoise.time_resample(tedges=master_time_post,update=True,remainder=False)
               
            
#            if tres_post > 0 or tres <= 0.5:
#                MolDenoise.time_resample(tedges=master_time_post,update=True,remainder=False)
        
        if settings['hsrl_rb_adjust']:
            print('Obtaining Rayleigh-Brillouin Correction')
            dnu = 20e6  # resolution
            nu_max = 10e9 # max frequency relative to line center
            nu = np.arange(-nu_max,nu_max,dnu)
            Ti2 = np.interp(nu,i2_data['freq']*1e9,i2_data['mol_scan'])  # molecular transmission
            
            Tc2 = np.interp(nu,i2_data['freq']*1e9,i2_data['combined_scan'])  # combined transmission
            
            [eta_i2,eta_c] = lp.RB_Efficiency([Ti2,Tc2],temp.profile.flatten(),pres.profile.flatten()*9.86923e-6,profs['molecular'].wavelength,nu=nu,norm=True,max_size=10000)
            
        #    beta_mol_norm = lp.RB_Spectrum(temp.profile.flatten(),pres.profile.flatten()*9.86923e-6,profs['molecular'].wavelength,nu=nu,norm=True)
        #    eta_i2 = np.sum(Ti2[:,np.newaxis]*beta_mol_norm,axis=0)
            eta_i2 = eta_i2.reshape(temp.profile.shape)
            profs['molecular'].multiply_piecewise(1.0/eta_i2)
        
            profs['molecular'].gain_scale(mol_gain,gain_var = (mol_gain*0.05)**2)
        
                
        #    eta_c = np.sum(Tc2[:,np.newaxis]*beta_mol_norm,axis=0)
            eta_c = eta_c.reshape(temp.profile.shape)
            profs['combined_hi'].multiply_piecewise(1.0/eta_c)
            profs['cross'].multiply_piecewise(1.0/eta_c)
            
            if settings['Denoise_Mol']:
                MolDenoise.multiply_piecewise(1.0/eta_i2)
                MolDenoise.gain_scale(mol_gain,gain_var = (mol_gain*0.05)**2)
                
            if settings['get_extinction'] and settings['as_altitude']:
                [eta_i2_ext] = lp.RB_Efficiency([Ti2],temp_ext.profile.flatten(),pres_ext.profile.flatten()*9.86923e-6,profs['molecular'].wavelength,nu=nu,norm=True,max_size=10000)
                eta_i2_ext = eta_i2_ext.reshape(temp_ext.profile.shape)
                mol_ext.multiply_piecewise(1.0/eta_i2_ext)
                mol_ext.range_correct()
                
        else:
            # Rescale molecular channel to match combined channel gain
            profs['molecular'].gain_scale(mol_gain,gain_var = (mol_gain*0.05)**2)
            if settings['Denoise_Mol']:
                MolDenoise.gain_scale(mol_gain)
        
        beta_a = lp.AerosolBackscatter(profs['molecular'],(profs['combined_hi']+profs['cross']),beta_m)
        if settings['Denoise_Mol']:
            beta_a_denoise = lp.AerosolBackscatter(MolDenoise,(profs['combined_hi']+profs['cross']),beta_m)
            beta_a_denoise.descript = 'Poisson total variation denoised calibrated measurement of Aerosol Backscatter Coefficient in m^-1 sr^-1'
            beta_a_denoise.label = 'Denoised Aerosol Backscatter Coefficient'
        
        if settings['get_extinction']:
            ext_sg_wid = settings['ext_sg_width']
            ext_sg_order = settings['ext_sg_order']
            ext_tres = settings['ext_tres']
            ext_zres = settings['ext_zres']
#            ext_sg_wid = 21
#            ext_sg_order = 4
#            ext_tres = 20  # extinction time resolution in seconds
#            ext_zres = 15   # extinction range resolution in meters
            mol_ext.conv(ext_tres/mol_ext.mean_dt,ext_zres/mol_ext.mean_dR)
            beta_m_ext.conv(ext_tres/beta_m_ext.mean_dt,ext_zres/beta_m_ext.mean_dR)
            OD = mol_ext/beta_m_ext
            OD.descript = 'Total optical depth from aircraft altitude'
            OD.label = 'Optical Depth'
            OD.profile_type = 'unitless'
            
            alpha_a = OD.copy()
            for ai in range(alpha_a.profile.shape[0]):
                alpha_a.profile[ai,:] = -2*gv.savitzky_golay(np.log(alpha_a.profile[ai,:].flatten()), ext_sg_wid, ext_sg_order, deriv=1)
            alpha_a = alpha_a/alpha_a.mean_dR # not sure this is the right scaling factor
            alpha_a = alpha_a - beta_m_ext*(8*np.pi/3)  # remove molecular extinction
            if settings['as_altitude']:
                alpha_a.range2alt(master_alt,air_data_post,telescope_direction=var_post['TelescopeDirection'])
                
    #        if tres_post > 0 or tres <= 0.5:
    #            alpha_a.time_resample(tedges=master_time_post,update=True,remainder=False)
            
            alpha_a.descript = 'Aerosol Extinction Coefficient'
            alpha_a.label = 'Aerosol Extinction Coefficient'
            alpha_a.profile_type = '$m^{-1}$'
        
        BSR = (profs['combined_hi']+profs['cross'])/profs['molecular']
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
        dVol.descript = 'Propensity of Volume to depolarize (d).  This is not identical to the depolarization ratio.  See Gimmestad: 10.1364/AO.47.003795 or Hayman and Thayer: 10.1364/JOSAA.29.000400'
        dVol.label = 'Volume Depolarization'
        dVol.profile_type = 'unitless'
        
        d_mol = 2*0.000365/(1+0.000365) # molecular depolarization
        
        #Particle Depolarization = dVol/(1.0-1.0/BSR) - d_mol/(BSR-1)
        dPart = (BSR*dVol-d_mol)/(BSR-1)
        dPart.descript = 'Propensity of Particles to depolarize (d).  This is not identical to the depolarization ratio.  See Gimmestad: 10.1364/AO.47.003795 or Hayman and Thayer: 10.1364/JOSAA.29.000400'
        dPart.label = 'Particle Depolarization'
        dPart.profile_type = 'unitless'
        
        
        
        if settings['Estimate_Mol_Gain']:
            # This segment estimates what the molecular gain should be 
            # based on a histogram minimum in BSR over the loaded data
            
            iUp = np.nonzero(var_post['TelescopeDirection']==1.0)[0]            
            lp.Estimate_Mol_Gain(BSR,iKeep=iUp,mol_gain=mol_gain_up,alt_lims=[2000,4000],label='Telescope Up',plot=True)
            
            
            iDown = np.nonzero(var_post['TelescopeDirection']==0.0)[0]
            lp.Estimate_Mol_Gain(BSR,iKeep=iDown,mol_gain=mol_gain_down,alt_lims=[2000,4000],label='Telescope Down',plot=True)
            
        
        
        
        # add a diagnostic for counts/backscatter coeff
        
        # add a diagnostic for diff overlap between lo and hi channels as a function
        # of count rate or backscatter coeff
        
        count_mask = profs['combined_hi'].profile < settings['count_mask_threshold']
        
        #dPartMask = dPart.profile_variance > 1.0
        dPartMask = dPart.profile_variance > settings['d_part_res_lim']**2
        #dPart.mask(dPartMask)
        dPart.mask(dPart.profile_variance > settings['d_part_res_lim']**2)
        dPart.mask(dPart.profile > 1.0)
        dPart.mask(dPart.profile < -0.1)
        
        try:
            proj_label = process_vars['proj_label']
        except KeyError:
            proj_label = ''
    #    proj_label = proj + ' ' + flight_label[usr_flt] + ', '
        
        beta_a.mask(np.isnan(beta_a.profile))
        BSR.mask(np.isnan(BSR.profile))
        dPart.mask(np.isnan(dPart.profile))
        dVol.mask(np.isnan(dVol.profile))
        
        if settings['Denoise_Mol']:
            beta_a_denoise.mask(np.isnan(beta_a_denoise.profile))
        
        
        if settings['count_mask_threshold'] > 0:
            beta_a.mask(count_mask)
            dPart.mask(count_mask)
            dVol.mask(count_mask)
            profs['combined_hi'].mask(count_mask)
            
            if settings['get_extinction']:
                alpha_a.mask(count_mask)
                alpha_a.mask(dPart.profile.mask)
        
        
        
        save_prof_list = [beta_a,dPart,dVol,BSR,beta_m]
        # add all channels to list of profilse to save
        for var in profs.keys():
            save_prof_list.extend(profs[var])
        if settings['get_extinction']:
            # add extinction to list of profilse to save
            save_prof_list.extend([alpha_a])
        if settings['Denoise_Mol']:
            # add denoised molecular observations to list of profilse to save
            save_prof_list.extend([MolDenoise])
        save_var1d_post = {'TelescopeDirection':{'description':'1-Lidar Pointing Up, 0-Lidar Pointing Down','units':'none'},
                           'polarization':{'description':'System Quarter Waveplate orientation','units':'degrees'}}
        save_air_post = {'THDG': {'description':'aircraft heading','units':'degrees'},
                         'TASX': {'description':'airspeed','units':'meters/second'},
                         'GGLAT': {'description':'latitude','units':'degrees'},
                         'PITCH': {'description':'aircraft pitch angle','units':'degrees'},
                         'GGALT': {'description':'altitude','units':'meters'},
                         'PSXC': {'description':'ambiant pressure','units':'hPa'},
                         'ROLL': {'description':'aircraft roll angle','units':'degrees'},
                         'GGLON': {'description':'longitude','units':'degrees'},
                         'ATX': {'description':'ambiant temperature', 'units':'C'}}
        
        if settings['save_data']:
            print('saving profiles')
            for ai in range(len(save_prof_list)):
                save_prof_list[ai].write2nc(save_data_file) #,name_override=True,tag=var_name)
                
            print('saving lidar status data')
            for var in save_var1d_post.keys():
                lp.write_var2nc(var_post[var],str(var),save_data_file,description=save_var1d_post[var]['description'],units=save_var1d_post[var]['units'])
            
            print('saving aircraft variables')    
            for var in save_air_post.keys():
                lp.write_var2nc(air_data_post[var],str(var),save_data_file,description=save_air_post[var]['description'],units=save_air_post[var]['units'])
        
        if settings['plot_2D']:
            if settings['plot_date']:
                t1d_plt = mdates.date2num([datetime.datetime.fromordinal(BSR.StartDate.toordinal()) \
                            + datetime.timedelta(seconds=sec) for sec in time_1d])   
            else:
                t1d_plt = time_1d/3600.0
            
            tlims = [(time_start-flight_date[usr_flt]).total_seconds()/3600.0,
                      (time_stop-flight_date[usr_flt]).total_seconds()/3600.0]
    #        tlims = [time_post[0]/3600.0, time_post[-1]/3600.0]
        #    rfig = lp.pcolor_profiles([BSR,dPart],scale=['log','linear'],climits=[[1,1e2],[0,0.7]],ylimits=[MinAlt*1e-3,MaxAlt*1e-3],title_add=proj_label,plot_date=plot_date)
            rfig = lp.pcolor_profiles([beta_a],scale=['log'],
                                      climits=[[1e-8,1e-3]],
                                      ylimits=[MinAlt*1e-3,MaxAlt*1e-3],
                                      tlimits=tlims,
                                      title_add=proj_label,
                                      plot_date=settings['plot_date'],
                                      t_axis_scale=settings['time_axis_scale'],
                                      h_axis_scale=settings['alt_axis_scale'],
                                      minor_ticks=5,major_ticks=1)
            if settings['as_altitude']:
                for ai in range(len(rfig[1])):
                    rfig[1][ai].plot(t1d_plt,air_data_t['GGALT']*1e-3,color='gray',linewidth=1.2)  # add aircraft altitude      
            if settings['save_plots']:
                plt.savefig(save_plots_path+'Aerosol_Backscatter_'+save_plots_base,dpi=300)
            
            rfig = lp.pcolor_profiles([dPart],scale=['linear'],
                                      climits=[[0,1.0]],
                                      ylimits=[MinAlt*1e-3,MaxAlt*1e-3],tlimits=tlims,
                                      title_add=proj_label,
                                      plot_date=settings['plot_date'],
                                      t_axis_scale=settings['time_axis_scale'],
                                      h_axis_scale=settings['alt_axis_scale'],
                                      minor_ticks=5,major_ticks=1)
            if settings['as_altitude']:
                for ai in range(len(rfig[1])):
                    rfig[1][ai].plot(t1d_plt,air_data_t['GGALT']*1e-3,color='gray',linewidth=1.2)  # add aircraft altitude      
            if settings['save_plots']:
                plt.savefig(save_plots_path+'Aerosol_Depolarization_'+save_plots_base,dpi=300)        
            
            rfig = lp.pcolor_profiles([profs['combined_hi']],scale=['log'],
                                      climits=[[1e-1,1e4]],
                                      ylimits=[MinAlt*1e-3,MaxAlt*1e-3],tlimits=tlims,
                                      title_add=proj_label,
                                      plot_date=settings['plot_date'],
                                      t_axis_scale=settings['time_axis_scale'],
                                      h_axis_scale=settings['alt_axis_scale'],
                                      minor_ticks=5,major_ticks=1)
            if settings['as_altitude']:
                for ai in range(len(rfig[1])):
                    rfig[1][ai].plot(t1d_plt,air_data_t['GGALT']*1e-3,color='gray',linewidth=1.2)  # add aircraft altitude
            if settings['save_plots']:
                plt.savefig(save_plots_path+'AttenuatedBackscatter_'+save_plots_base,dpi=300)
            
            rfig = lp.pcolor_profiles([dVol],scale=['linear'],
                                      climits=[[0,1.0]],
                                      ylimits=[MinAlt*1e-3,MaxAlt*1e-3],tlimits=tlims,
                                      title_add=proj_label,
                                      plot_date=settings['plot_date'],
                                      t_axis_scale=settings['time_axis_scale'],
                                      h_axis_scale=settings['alt_axis_scale'],
                                      minor_ticks=5,major_ticks=1)
            if settings['as_altitude']:
                for ai in range(len(rfig[1])):
                    rfig[1][ai].plot(t1d_plt,air_data_t['GGALT']*1e-3,color='gray',linewidth=1.2)  # add aircraft altitude
            if settings['save_plots']:
                plt.savefig(save_plots_path+'Volume_Depolarization_'+save_plots_base,dpi=300)
                    
            
            if settings['get_extinction']:
                rfig = lp.pcolor_profiles([alpha_a],scale=['log'],
                                          climits=[[1e-5,1e-2]],
                                          ylimits=[MinAlt*1e-3,MaxAlt*1e-3],tlimits=tlims,
                                          title_add=proj_label,
                                          plot_date=settings['plot_date'],
                                          t_axis_scale=settings['time_axis_scale'],
                                          h_axis_scale=settings['alt_axis_scale'],
                                          minor_ticks=5,major_ticks=1)
                if settings['as_altitude']:
                    for ai in range(len(rfig[1])):
                        rfig[1][ai].plot(t1d_plt,air_data_t['GGALT']*1e-3,color='gray',linewidth=1.2)  # add aircraft altitude
                if settings['save_plots']:
                    plt.savefig(save_plots_path+'Extinction_'+save_plots_base,dpi=300)
                    
            if settings['Denoise_Mol']:
                rfig = lp.pcolor_profiles([beta_a_denoised],scale=['log'],
                                      climits=[[1e-8,1e-3]],
                                      ylimits=[MinAlt*1e-3,MaxAlt*1e-3],
                                      tlimits=tlims,
                                      title_add=proj_label,
                                      plot_date=settings['plot_date'],
                                      t_axis_scale=settings['time_axis_scale'],
                                      h_axis_scale=settings['alt_axis_scale'],
                                      minor_ticks=5,major_ticks=1)
                if settings['as_altitude']:
                    for ai in range(len(rfig[1])):
                        rfig[1][ai].plot(t1d_plt,air_data_t['GGALT']*1e-3,color='gray',linewidth=1.2)  # add aircraft altitude      
                if settings['save_plots']:
                    plt.savefig(save_plots_path+'Denoised_Aerosol_Backscatter_'+save_plots_base,dpi=300)
                
                
                
            #lp.plotprofiles(profs)
            #dPart.mask(dPartMask)
            #lp.pcolor_profiles([BSR,dVol],scale=['log','linear'],climits=[[1,5e2],[0,1.0]])
            #lp.pcolor_profiles([dVol],scale=['linear'],climits=[[0,1]])
        if settings['show_plots']:
            plt.show()
            
        return save_prof_list
    else:
        # no data was loaded
        return 0