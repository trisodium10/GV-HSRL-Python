# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 11:38:22 2018

@author: mhayman
"""
import os
import sys
import datetime
import numpy as np

import matplotlib
#matplotlib.use('Agg')  # disable interactive displays (allows running without xwindow)

import matplotlib.pyplot as plt



"""
Process a full flight in hour increments
"""
process_start_time = datetime.datetime.now()
process_vars = {}
process_vars['proj'] = 'SOCRATES'

#process_vars['flt'] = process_vars['proj']+'rf07'
##process_vars['flight_time_start'] = datetime.timedelta(hours=6,minutes=25)
#process_vars['flight_time_start'] = datetime.timedelta(hours=2,minutes=00)
##process_vars['flight_time_start'] = datetime.timedelta(hours=25,minutes=5)
#process_vars['flight_time_stop'] = process_vars['flight_time_start'] + datetime.timedelta(hours=0,minutes=25)

#process_vars['flt'] = process_vars['proj']+'rf12'
#process_vars['flight_time_start'] = datetime.timedelta(hours=25,minutes=0) # 24:55
###process_vars['flight_time_start'] = datetime.timedelta(hours=25,minutes=5)
#process_vars['flight_time_stop'] = datetime.timedelta(hours=26,minutes=0) # 25:15
###process_vars['flight_time_stop'] = process_vars['flight_time_start'] + datetime.timedelta(hours=0,minutes=30)

process_vars['flt'] = process_vars['proj']+'rf02'
process_vars['flight_time_start'] = datetime.timedelta(hours=5,minutes=0) # 24:55
process_vars['flight_time_stop'] = datetime.timedelta(hours=6,minutes=0) # 25:15

# size of each processing step
time_increment = datetime.timedelta(hours=1,minutes=0)
# size of a processesed data set
time_duration = datetime.timedelta(hours=1,minutes=0)
#time_duration = time_increment

settings = {
    'full_flight':False, # process the entire flight
    'tres':0.5,  # resolution in time in seconds (0.5 sec) before altitude correction
    'tres_post':2.0, # resolution after altitude correction (in seconds) -  set to zero to not use
    'zres':7.5,  # altitude resolution in meters (7.5 m minimum)

    'mol_smooth':False, # smooth molecular profile
    
    #mol_gain = 1.133915#1.0728915  # gain adjustment to molecular channel
    
    # index for where to treat the profile as background only
    'BGIndex': -100, # negative number provides an index from the end of the array
    'platform':'airborne', # 'ground' or 'airborne'.  If 'airborne' it needs an aircraft netcdf.
    'MaxAlt':14e3,
    'MinAlt':-1e3,
    
    'RemoveCals':True,  # don't include instances where the I2 cell is removed
                        # scan files are not included in the file search so they
                        # are removed anyway
    
    'Remove_Off_Data':True, # remove instances where the lidar does not appear
                            # to be running
    
    'get_extinction':True, # retrieve extinction estimate
    
    'diff_geo_correct':True,  # apply differential overlap correction
    'deadtime_correct':True,  # perform deadtime correction
    
    'load_reanalysis':False, # load T and P reanalysis from NCEP/NCAR Model
    
    'plot_2D':True,   # pcolor plot the BSR and depolarization profiles
    'show_plots':True, # show plots in a matplotlib window
    'plot_date':True,  # plot results in date time format.  Otherwise plots as hour floats
    
    'save_plots':False, # save the plot data
    
    'save_data':False, # save data as netcdf
    
    'save_flight_folder':True, # save data/plots in folders according to flight name
    
    'time_axis_scale':5.0,  # scale for horizontal axis on pcolor plots    
    'count_mask_threshold':2.0,  # count mask threshold (combined_hi).  If set to zero, no mask applied 
    
    'Estimate_Mol_Gain':True, # use statistics on BSR to estimate the molecular gain
    'save_mol_gain_plot':False,  # save the plots used to estimate molecular gain
    
    'hsrl_rb_adjust':True, # adjust for Rayleigh Brillouin Spectrum
    
    'Denoise_Mol':False, # run PTV denoising on molecular channel
    'denoise_accel':True, # run accelerated denoising 
        
    'time_denoise':False,  # run horizontal (time) denoising on profiles
    'time_denoise_accel':True, # run accelerated denoising (reduced scan region)
    'time_denoise_debug_plots':False, # plot denoising results (use only for debugging.  Slows processing down.)
    'time_denoise_eps':1e-5,  # eps float precision for optimization.  Smaller numbers are slower but provide more accurate denoising 
    'time_denoise_verbose':False,  # output optimizor status at each step
    'time_denoise_max_range':3e3,  # maximum range to which denoising occurs
    
    'Airspeed_Threshold':25, # threshold for determining start and end of the flight (in m/s)
    
    'loadQWP':'fixed',  # load 'fixed','rotating', or 'all' QWP data
    
    'as_altitude':True, # process in altitude centered format or range centered format

    'aircraft_time_shift':0.859822  # shift in aircraft time needed to align to HSRL time 0.75 - 3.0
    }
    
    


PathFile = os.path.abspath(__file__+'/../')+'/gv_hsrl_socrates_paths.py'

print('Paths stored in: ' +PathFile)

# load path data for this computer
exec(open(PathFile).read())

# add the path to GVHSRLlib manually
library_path = os.path.abspath(paths['software_path']+'/processors/')
print(library_path)
if library_path not in sys.path:
    sys.path.append(library_path)

import Airborne_GVHSRL_DataProcessor as dp
import Airborne_GVHSRL_DataSelection as ds

#Processor = software_path + 'processors/Airborne_GVHSRL_DataProcessor.py'
#DataSelector = software_path + 'processors/Airborne_GVHSRL_DataSelection.py'

time_start0,time_stop0,settings,paths,process_vars = \
    ds.SelectAirborneData(settings=settings,paths=paths,process_vars=process_vars)



day_start = datetime.datetime(year=time_start0.year,month=time_start0.month,day=time_start0.day)
t0 = time_start0-day_start

print(time_start0)
print(time_stop0)

prof_list = dp.ProcessAirborneDataChunk(time_start0,time_stop0,
                             settings=settings,paths=paths,process_vars=process_vars,date_reference=day_start)


#time_start = time_start0
#time_stop = day_start + datetime.timedelta(seconds=np.floor(t0.total_seconds()/time_increment.total_seconds())*time_increment.total_seconds())+time_duration
#run_loop = False
#while run_loop:
##    # check if we are going to overrrun the data
##    if time_stop > time_stop0:
##        time_stop = time_stop0
#        
##    g['time_start'] = time_start
##    g['time_stop'] = time_stop
##    
##    exec(open(Processor).read())
#        
#    # only save data every hour increment
##    if time_start.minute == 0:
##        settings['save_data'] = True
#    if time_stop >= time_stop0:
#        run_loop = False
##    else:
##        settings['save_data'] = False
#        
#    
#    dp.ProcessAirborneDataChunk(time_start,time_stop,
#                             settings=settings,paths=paths,process_vars=process_vars,date_reference=day_start)
#
#    run_loop = False
#    
#
#    time_start = time_stop
#    time_stop = time_stop+time_increment
#    #    time_start = time_stop #+datetime.timedelta(seconds=0.1)
process_stop_time = datetime.datetime.now()
plt.show()

print('processing began ' + process_start_time.strftime('%Y-%b-%d %H:%M'))
print('processing completed ' + process_stop_time.strftime('%Y-%b-%d %H:%M'))
print('processing duration %f hours' %((process_stop_time-process_start_time).total_seconds()/3600.0))
