# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 11:38:22 2018

@author: mhayman
"""
import os
import datetime
import numpy as np



"""
Process a full flight in hour increments
"""


proj = 'SOCRATES'
#flt = proj+'tf02'
#flight_time_start = datetime.timedelta(hours=18,minutes=00)
#flight_time_stop = datetime.timedelta(hours=19,minutes=0)
#flight_time_stop = flight_time_start + datetime.timedelta(hours=0,minutes=30)

# size of each processing step
time_increment = datetime.timedelta(hours=1,minutes=00)

settings = {
    'tres':0.5,  # resolution in time in seconds (0.5 sec) before altitude correction
    'tres_post':10.0, # resolution after altitude correction (in seconds) -  set to zero to not use
    'zres':7.5,  # altitude resolution in meters (7.5 m minimum)
    
    #mol_gain = 1.133915#1.0728915  # gain adjustment to molecular channel
    
    # index for where to treat the profile as background only
    'BGIndex': -100, # negative number provides an index from the end of the array
    'platform':'airborne', # 'ground' or 'airborne'.  If 'airborne' it needs an aircraft netcdf.
    'MaxAlt':14e3,
    'MinAlt':-30,
    
    'RemoveCals':True,  # don't include instances where the I2 cell is removed
                        # scan files are not included in the file search so they
                        # are removed anyway
    
    'Remove_Off_Data':True, # remove instances where the lidar does not appear
                            # to be running
    
    'diff_geo_correct':True,  # apply differential overlap correction
    
    'load_reanalysis':False, # load T and P reanalysis from NCEP/NCAR Model
    
    'plot_2D':True,   # pcolor plot the BSR and depolarization profiles
    'show_plots':False, # show plots in a matplotlib window
    'plot_date':True,  # plot results in date time format.  Otherwise plots as hour floats
    
    'save_plots':True, # save the plot data
    
    'save_data':False, # save data as netcdf
    
    'time_axis_scale':5.0,  # scale for horizontal axis on pcolor plots    
    'count_mask_threshold':2.0,  # count mask threshold (combined_hi).  If set to zero, no mask applied 
    
    'Estimate_Mol_Gain':False, # use statistics on BSR to estimate the molecular gain
    
    'hsrl_rb_adjust':True, # adjust for Rayleigh Brillouin Spectrum
    
    'Denoise_Mol':False, # run PTV denoising on molecular channel
    
    
    'Airspeed_Threshold':15, # threshold for determining start and end of the flight (in m/s)
    
    'loadQWP':'fixed',  # load 'fixed','rotating', or 'all' QWP data
    
    'as_altitude':True # process in altitude centered format or range centered format
    }
    
    
##basepath = '/scr/eldora1/HSRL_data/'  # old path - still works with link from HSRL_data to /hsrl/raw/
#basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path
##basepath = '/Users/mhayman/Documents/HSRL/GVHSRL_data/'  # local computer data path
        
#aircraft_basepath = {
#    'CSET':'/scr/raf_data/CSET/24_Mar_17_BU/',
#    'SOCRATES':'/scr/raf_data/SOCRATES/' #SOCRATEStf01.nc       
#    } 
    
#basepath = '/scr/rain1/rsfdata/projects/socrates/hsrl/raw/'
#
#save_data_path = os.environ['HOME'] + '/Documents/HSRL/GVHSRL/data/'
#save_plots_path = os.environ['HOME'] + '/Documents/HSRL/GVHSRL/plots/'

#Processor = os.environ['HOME']+'/Documents/Python/Lidar/GV-HSRL-Python/processors/Airborne_Processor_GVHSRL.py'

PathFile = os.path.abspath(__file__+'/../')+'/gv_hsrl_socrates_paths.py'



# load path data for this computer
exec(open(PathFile).read())



Processor = software_path + 'processors/Airborne_GVHSRL_DataProcessor.py'
DataSelector = software_path + 'processors/Airborne_GVHSRL_DataSelection.py'
fullpath = os.path.abspath(Processor)
#g = globals().copy()
#g['__file__'] = fullpath

# obtain the time references for processing the flight
exec(open(DataSelector).read())

time_start0 = time_start
time_stop0 = time_stop
day_start = datetime.datetime(year=time_start0.year,month=time_start0.month,day=time_start0.day)
t0 = time_start0-day_start

time_start = time_start0
time_stop = day_start + datetime.timedelta(seconds=np.floor(t0.total_seconds()/time_increment.total_seconds())*time_increment.total_seconds())+time_increment

while time_start <= time_stop0:
    # check if we are going to overrrun the data
    if time_stop > time_stop0:
        time_stop = time_stop0
        
    exec(open(Processor).read())
    
    time_start = time_stop+datetime.timedelta(seconds=0.1)
    time_stop = time_stop+time_increment
    
    