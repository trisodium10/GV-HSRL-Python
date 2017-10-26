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



# input and raw_input are the same for python 3 and 2 respectively
# this makes it so input always accepts a string
try:
    input=raw_input
except NameError:
    pass

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

tres = 1*60.0  # resolution in time points (0.5 sec)
zres = 10.0  # resolution in altitude points (7.5 m)

# index for where to treat the profile as background only
BGIndex = -100; # negative number provides an index from the end of the array
platform = 'ground' # 'ground' or 'airborne'.  If 'airborne' it needs an aircraft netcdf.
MaxAlt = 10e3


pol_xtalk = 0.015

#kB = 1.3806504e-23;
#c = 3e8


FilterI2 = True  # only include data where the I2 cell is removed


var_1d_list = ['total_energy','RemoveLongI2Cell'\
    ,'TelescopeDirection','TelescopeLocked','polarization']  # 'DATA_shot_count'
var_1d_data = dict(zip(var_1d_list,[np.array([])]*len(var_1d_list)))
timeD = np.array([])

#var_2d_list = ['cross','molecular','combined_hi','combined_lo']
#var_2d_list

#basepath = '/scr/eldora1/HSRL_data/'  # old path - still works with link from HSRL_data to /hsrl/raw/
basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path

#if Day < 10:
#    DayStr = '0' + str(Day)
#else:
#    DayStr = str(Day)
#    
#if Month < 10:
#    MonthStr = '0' + str(Month)
#else:
#    MonthStr = str(Month)
#
#YearStr = str(Year)

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
    


firstFile = True
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
        
        for ivar in range(len(var_1d_data.keys())):
            if any(var_1d_data.keys()[ivar] in s for s in f.variables):
                var_data = lp.ncvar(f,var_1d_data.keys()[ivar])
                var_1d_data[var_1d_data.keys()[ivar]] = np.concatenate((var_1d_data[var_1d_data.keys()[ivar]],var_data)) 
        
        
        
#        #TransEnergy = f.variables['total_energy'].data.copy();      # Transmitted pulse energy
#        TransEnergy0 = lp.ncvar(f,'total_energy')
#        
#        cell_status = lp.ncvar(f,'RemoveLongI2Cell')        
#        
        # System for time array still needs work
        timeD0 = np.array(f.variables['DATA_time'][:])
        time_dt0 = [datetime.datetime(*x) for x in timeD0[:,:7]]
        time_sec0 = np.array([(x-time_dt0[0]).total_seconds() for x in time_dt0])
        
        if len(timeD)>0:
            timeD = np.vstack((timeD,np.array(f.variables['DATA_time'][:])));          # time array [year,month,day,hour,minute,second,msec,usec]
        else:
            timeD = np.array(f.variables['DATA_time'][:])
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
        
        print(f.variables['cross'][:].shape())     
        
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
time_dt = [datetime.datetime(*x) for x in timeD[:,:7]]
time_sec = np.array([(x-time_dt[0]).total_seconds() for x in time_dt])



#        # load profile data
#        if firstFile:
#            TransEnergy = TransEnergy0.copy()
#            timeDsec = timeDsec0.copy()
#            shot_count = shot_count0.copy()
#            telescope_direction = telescope_direction0.copy()
#            telescope_locked = telescope_locked0.copy()            
#            
#            CrossPol = lp.LidarProfile(cross_data,timeDsec0,label='Cross Polarization Channel',descript = 'Cross Polarization\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')
#            RemCross = CrossPol.time_resample(delta_t=tres,update=True,remainder=True)
#            
#            Molecular = lp.LidarProfile(mol_data,timeDsec0,label='Molecular Backscatter Channel',descript = 'Parallel Polarization\nMolecular Backscatter Returns',bin0=47,lidar='GV-HSRL')
#            RemMol = Molecular.time_resample(delta_t=tres,update=True,remainder=True)
#            
#            CombHi = lp.LidarProfile(hi_data,timeDsec0,label='High Gain Total Backscatter Channel',descript = 'Parallel Polarization\nHigh Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')
#            RemCombHi = CombHi.time_resample(delta_t=tres,update=True,remainder=True)
#            
#            CombLo = lp.LidarProfile(lo_data,timeDsec0,label='Low Gain Total Backscatter Channel',descript = 'Parallel Polarization\nLow Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')            
#            RemCombLo = CombLo.time_resample(delta_t=tres,update=True,remainder=True)
#            
#            
##            Molecular = lp.LidarProfile(mol_data.T,dt*np.arange(mol_data.shape[1])+Hour*3600,label='Molecular Backscatter Channel',descript = 'Unpolarization\nMolecular Backscatter Returns',bin0=0,lidar='DLB-HSRL')
##            RemMol = Molecular.time_resample(i=tres,update=True,remainder=True)
##            
##            CombHi = lp.LidarProfile(hi_data.T,dt*np.arange(hi_data.shape[1])+Hour*3600,label='Total Backscatter Channel',descript = 'Unpolarization\nHigh Gain\nCombined Aerosol and Molecular Returns',bin0=0,lidar='DLB-HSRL')
##            RemCom = CombHi.time_resample(i=tres,update=True,remainder=True)
#            
#            firstFile = False
#            
#            
#        else:
#            TransEnergy = np.concatenate((TransEnergy,TransEnergy0))
#            timeDsec = np.concatenate((timeDsec,timeDsec0))
#            shot_count = np.concatenate((shot_count,shot_count0))
#            telescope_direction = np.concatenate((telescope_direction,telescope_direction0))
#            telescope_locked = np.concatenate((telescope_locked,telescope_locked0))
#            
#            if np.size(RemMol.time) > 0:
#                MolTmp = lp.LidarProfile(mol_data,timeDsec0,label='Molecular Backscatter Channel',descript = 'Parallel Polarization\nMolecular Backscatter Returns',bin0=47,lidar='GV-HSRL')
#                MolTmp.cat_time(RemMol)
#                RemMol = MolTmp.time_resample(delta_t=tres,update=True,remainder=True)
#                Molecular.cat_time(MolTmp,front=False)
#                
#                CrossTmp = lp.LidarProfile(cross_data,timeDsec0,label='Cross Polarization Channel',descript = 'Cross Polarization\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')
#                CrossTmp.cat_time(RemCross)
#                RemCross = CrossTmp.time_resample(delta_t=tres,update=True,remainder=True)
#                CrossPol.cat_time(CrossTmp,front=False)
#                
#                CombHiTmp = lp.LidarProfile(hi_data,timeDsec0,label='High Gain Total Backscatter Channel',descript = 'Parallel Polarization\nHigh Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')
#                CombHiTmp.cat_time(RemCombHi)
#                RemCombHi = CombHiTmp.time_resample(delta_t=tres,update=True,remainder=True)
#                CombHi.cat_time(CombHiTmp,front=False)
#                
#                CombLoTmp = lp.LidarProfile(lo_data,timeDsec0,label='Low Gain Total Backscatter Channel',descript = 'Parallel Polarization\nLow Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')            
#                CombLoTmp.cat_time(RemCombHi)
#                RemCombLo = CombLoTmp.time_resample(delta_t=tres,update=True,remainder=True)
#                CombLo.cat_time(CombLoTmp,front=False)
#            else:
#                MolTmp = lp.LidarProfile(mol_data,timeDsec0,label='Molecular Backscatter Channel',descript = 'Parallel Polarization\nMolecular Backscatter Returns',bin0=47,lidar='GV-HSRL')
#                RemMol = MolTmp.time_resample(i=tres,update=True,remainder=True)
#                Molecular.cat_time(MolTmp,front=False)
#                
#                CrossTmp = lp.LidarProfile(cross_data,timeDsec0,label='Cross Polarization Channel',descript = 'Cross Polarization\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')
#                RemCross = CrossTmp.time_resample(delta_t=tres,update=True,remainder=True)
#                CrossPol.cat_time(CrossTmp,front=False)
#                
#                CombHiTmp = lp.LidarProfile(hi_data,timeDsec0,label='High Gain Total Backscatter Channel',descript = 'Parallel Polarization\nHigh Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')
#                RemCombHi = CombHiTmp.time_resample(delta_t=tres,update=True,remainder=True)
#                CombHi.cat_time(CombHiTmp,front=False)
#                
#                CombLoTmp = lp.LidarProfile(lo_data,timeDsec0,label='Low Gain Total Backscatter Channel',descript = 'Parallel Polarization\nLow Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')            
#                RemCombLo = CombLoTmp.time_resample(delta_t=tres,update=True,remainder=True)
#                CombLo.cat_time(CombLoTmp,front=False)