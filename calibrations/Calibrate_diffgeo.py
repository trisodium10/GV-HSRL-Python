# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 09:25:01 2017

@author: mhayman
"""
import numpy as np
import matplotlib.pyplot as plt
#from scipy.io import netcdf
import netCDF4 as nc4
import LidarProfileFunctions as lp
import scipy.interpolate
import datetime
import glob



#def match_data_times(master_time,profiles,var_1d_data,time_sec):
#    for iprof in range(len(profiles)):
#        profiles[iprof].time_resample(tedges=master_time,update=True,remainder=False)
#    
#    itime = np.digitize(time_sec,master_time)
    


# input and raw_input are the same for python 3 and 2 respectively
# this makes it so input always accepts a string
try:
    input=raw_input
except NameError:
    pass


year_in = 2017
month_in = 10
day_in = 25
start_hr = 15
stop_hr = 6

print('Default Date:')
print('(M/D/Y) %d/%d/%d, starting %.1f UTC for %.1f h'%(month_in,day_in,year_in,start_hr,stop_hr))
if input('Run this default date? [y/n]') != 'y':
    print("Enter search range for diff_geo calibration:")
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

tres = 1*60.0  # resolution in time in seconds (0.5 sec)
zres = 10.0  # resolution in altitude points (7.5 m)

# index for where to treat the profile as background only
BGIndex = -100; # negative number provides an index from the end of the array
platform = 'ground' # 'ground' or 'airborne'.  If 'airborne' it needs an aircraft netcdf.
MaxAlt = 10e3


pol_xtalk = 0.015

#kB = 1.3806504e-23;
#c = 3e8


FilterI2 = False  # only include data where the I2 cell is removed


var_1d_list = ['total_energy','RemoveLongI2Cell'\
    ,'TelescopeDirection','TelescopeLocked','polarization']  # 'DATA_shot_count'
var_1d_data = dict(zip(var_1d_list,[np.array([])]*len(var_1d_list)))
timeD = np.array([])

var_2d_list = ['cross','molecular','combined_hi','combined_lo']
var_2d_data = dict(zip(var_2d_list,[np.array([])]*len(var_2d_list)))


#basepath = '/scr/eldora1/HSRL_data/'  # old path - still works with link from HSRL_data to /hsrl/raw/
basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path


#FilePath0 = basepath + YearStr + '/' + MonthStr + '/' + DayStr + '/raw/'
day_iter =  cal_start.date()
SubFiles = []
FileDate = []
while day_iter <= cal_stop.date():
    FilePath0 = basepath+day_iter.strftime('%Y/%m/%d/raw/')
    SubList0 = glob.glob(FilePath0+'*data*.nc')
    SubFiles = SubFiles+SubList0                  # make a list of files falling on the dates corresponding to the search period
    FileDate = FileDate+len(SubList0)*[day_iter]  # collect dates associated with each file (it's just easier)
    day_iter = day_iter + datetime.timedelta(days=1)
    


#firstFile = True
for idir in range(len(SubFiles)):
    # get the file start time from the file name
    iTime = np.char.find(SubFiles[idir],'T')
    Hour = np.double(SubFiles[idir][iTime+1:iTime+3])
    File_Time = datetime.datetime(*FileDate[idir].timetuple()[:4])+datetime.timedelta(hours=Hour)
    # check if this file is in the search time period
    if File_Time >= cal_start and File_Time <= cal_stop:
#        if Hour >= np.floor(HourLim[0]) and Hour <= HourLim[1]:
        # Open the netcdf file and load the relevant variables
#        f = netcdf.netcdf_file(SubFiles[idir], 'r')
        
        f = nc4.Dataset(SubFiles[idir],'r')
        
        for var in var_1d_data.keys():
            if any(var in s for s in f.variables):
                var_data = lp.ncvar(f,var)
#                print(var+' length: %d'var_data.size)
                var_1d_data[var] = np.concatenate((var_1d_data[var],var_data)) 
            
#        for ivar in range(len(var_1d_data.keys())):
#            if any(list(var_1d_data.keys())[ivar] in s for s in f.variables):
#                var_data = lp.ncvar(f,list(var_1d_data.keys())[ivar])
#                var_1d_data[list(var_1d_data.keys())[ivar]] = np.concatenate((var_1d_data[list(var_1d_data.keys())[ivar]],var_data)) 
        
        
        
#        #TransEnergy = f.variables['total_energy'].data.copy();      # Transmitted pulse energy
#        TransEnergy0 = lp.ncvar(f,'total_energy')
#        
#        cell_status = lp.ncvar(f,'RemoveLongI2Cell')        
#        
        # System for time array still needs work
        timeD0 = np.array(f.variables['DATA_time'][:]).astype(np.int)
        timeD0[:,6] = timeD0[:,6]*1e3+timeD0[:,7]
        time_dt0 = [datetime.datetime(*x) for x in timeD0[:,:7]]
        time_sec0 = np.array([(x-time_dt0[0]).total_seconds() for x in time_dt0])
        
        if len(timeD)>0:
            timeD = np.vstack((timeD,np.array(f.variables['DATA_time'][:].astype(np.int))));          # time array [year,month,day,hour,minute,second,msec,usec]
        else:
            timeD = np.array(f.variables['DATA_time'][:]).astype(np.int)
#        timeDsec0 = timeD[:,3]*3600.0+timeD[:,4]*60.0+timeD[:,5]+timeD[:,6]*1e-3;    # convert time to seconds starting at the hour
#        
#        shot_count0 = f.variables['DATA_shot_count'].data.copy()     # number of laser shots in the time bin
#        
#        telescope_direction0 = f.variables['TelescopeDirection'].data.copy()  # indicates if telescope is pointed up or down (1 = up, -1? = down)
#        telescope_locked0 = f.variables['TelescopeLocked'].data.copy()  # indicates if the telescope is locked - indicates possible tilted operation
#
#        # Check if polarization (QWP Rotation) is in this netcdf.
        if any('polarization' in s for s in f.variables):
            QWP = f.variables['polarization'][:].copy()               # QWP rotation angle
            # Calculate the QWP rotation rate to determine its status
            meanQWPrate = np.median(np.diff(QWP)/np.diff(time_sec0)*180/np.pi)
            meanQWPvalue = np.median(QWP*180/np.pi)
            if meanQWPrate >= 10:
                QWP_Status = 'rotating'
            else:
                QWP_Status = 'fixed'
            print( 'QWP Rotation Rate: %f deg/sec' %meanQWPrate)
            print( 'Average QWP Position: %f deg' %meanQWPvalue)
        else:
            QWP_Status = 'fixed'
        
        
        print( 'Processing %d UT' %Hour)
        print( 'Profile Integration Time: %f seconds' %np.median(np.diff(time_sec0)))
        print( 'Processing QWP as %s' %QWP_Status)
        
        for var in var_2d_data.keys():
            if len(var_2d_data[var]) == 0:
                var_2d_data[var] = f.variables[var][:].copy()
            else:
                var_2d_data[var] = np.vstack((var_2d_data[var],f.variables[var][:]))
                
        
#        print(f.variables['cross'][:].copy().shape)     
        
#        if len(cross_data) != 0:
#            cross_data = np.hstack()
#            
#        cross_data = lp.profile_netcdf(f.variables['cross'])
#        mol_data = lp.profile_netcdf(f.variables['molecular'])
#        hi_data = lp.profile_netcdf(f.variables['combined_hi'])
#        lo_data = lp.profile_netcdf(f.variables['combined_lo'])
        
##        print('%d,%d'(hi_data.data.shape[0],hi_data.data.shape[1]))
        f.close()
#timeD = timeD[:,:8]
timeD[:,6] = timeD[:,6]*1e3+timeD[:,7]
time_dt = [datetime.datetime(*x) for x in timeD[:,:7]]
time_sec = np.array([(x-time_dt[0]).total_seconds() for x in time_dt])

time_dt = np.array(time_dt)

plt.figure(); 
plt.plot(time_sec/3600,var_1d_data['RemoveLongI2Cell'])
plt.grid(b=True)
plt.xlabel('time [h-UTC]')
plt.ylabel('Value')
plt.title('RemoveLongI2Cell')

master_time = np.arange(time_sec[0]-tres/2,time_sec[-1]+tres/2,tres)

#match_data_times(master_time,profiles,var_1d_data,time_sec)
if 'cross' in var_2d_data.keys():
    CrossPol = lp.LidarProfile(var_2d_data['cross'],time_sec,label='Cross Polarization Channel',descript = 'Cross Polarization\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL',StartDate=cal_start.date())
    CrossPol.time_resample(tedges=master_time,update=True,remainder=False)
    CrossPol.bg_subtract(BGIndex)

if 'molecular' in var_2d_data.keys():
    Molecular = lp.LidarProfile(var_2d_data['molecular'],time_sec,label='Molecular Backscatter Channel',descript = 'Parallel Polarization\nMolecular Backscatter Returns',bin0=47,lidar='GV-HSRL',StartDate=cal_start.date())
    Molecular.time_resample(tedges=master_time,update=True,remainder=False)
    Molecular.bg_subtract(BGIndex)

if 'combined_hi' in var_2d_data.keys():
    CombHi = lp.LidarProfile(var_2d_data['combined_hi'],time_sec,label='High Gain Total Backscatter Channel',descript = 'Parallel Polarization\nHigh Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL',StartDate=cal_start.date())
    CombHi.time_resample(tedges=master_time,update=True,remainder=False)
    CombHi.bg_subtract(BGIndex)

if 'combined_lo' in var_2d_data.keys():       
    CombLo = lp.LidarProfile(var_2d_data['combined_lo'],time_sec,label='Low Gain Total Backscatter Channel',descript = 'Parallel Polarization\nLow Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL',StartDate=cal_start.date())            
    CombLo.time_resample(tedges=master_time,update=True,remainder=False)
    CombLo.bg_subtract(BGIndex)

if FilterI2:
    i2_size = var_1d_data['RemoveLongI2Cell'].size  # size of array.  Don't apply to any arrays that don't have matching time dimensions
    i2_rem = np.nonzero(var_1d_data['RemoveLongI2Cell'] < 50)[0]  # data points where the I2 cell is removed
    
    for var in var_1d_data.keys():
            if var_1d_data[var].size == i2_size:
                var_1d_data[var] = var_1d_data[var][i2_rem]
            else:
                print(var+' does not have time dimension that matches RemoveLongI2Cell')
#    for var in var_2d_data.keys():
#            if var_2d_data[var].shape[0] == i2_size:
#                var_2d_data[var] = var_2d_data[var][i2_rem,:]
#            else:
#                print(var+' does not have time dimension that matches RemoveLongI2Cell')
    
    time_dt = time_dt[i2_rem]
    time_sec = time_sec[i2_rem]
    
fig_data = lp.pcolor_profiles([CombHi,Molecular],ylimits=[0,12])  
fig_data[1][1].plot(time_sec/3600,var_1d_data['RemoveLongI2Cell']/15,'--')      

if all(var_1d_data['RemoveLongI2Cell'] > 50):
    print('No intervals found with I2 cell removed')
    RunCal = False

else:
    RunCal = True
    # find times where the i2 cell was inserted/removed
    ical = np.nonzero(np.diff(var_1d_data['RemoveLongI2Cell'])!=0)[0]
    if var_1d_data['RemoveLongI2Cell'][0] > 50:
        i0 = ical[::2]  # start indices
        i1 = ical[1::2] # stop indices
    else:
        i0 = np.concatenate((np.zeros(1),ical[1::2]))  # start indices
        i1 = ical[::2] # stop indices
        
    # ask user to select the interval to use
    print('Found calibration intervals:')
    for ai in range(len(i0)):
        if ai < i1.size:
            print('%d.)  %.2f - %.2f UTC'%(ai,time_sec[i0[ai]]/3600,time_sec[i1[ai]]/3600))
        else:
            print('%d.)  %.2f - End of file UTC'%(ai,time_sec[i0[ai]]/3600))
    cal_index = np.int(input('Select Interval (invalid number quits cal)'))
    

# check for power stability
# add 5 min buffers on either side
# trim profiles to the selected interval

# perform diff_geo

