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

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
#    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
    
def fit_high_range(profile,profile_var,i_min,i_max,order=1):
    """
    returns a profile with a polynomial fit to data between specified
    integers.
    profile - the profile to fit
    profile_var - variance of the supplied profile
    i_min - start of data used for the fit
    i_max - end of data used for the fit
    order - order of the polynomial fit (defaults to linear)
    """    
    
    xfit = np.arange(profile[i_min:i_max].size)
    yfit = profile[i_min:i_max]
    wfit = 1.0/np.sqrt(profile_var[0,i_min:i_max].flatten())
    wfit[0:5] = 10*np.max(wfit)
    pfit = np.polyfit(xfit,yfit,order,w=wfit)
    xprof = np.arange(hi_smooth[i_min:].size)
    profile_out = profile.copy()
    profile_out[i_const:] = np.polyval(pfit,xprof)    
    return profile_out

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
start_hr = 17
stop_hr = 4

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
time_dt0 = datetime.datetime(time_dt[0].year,time_dt[0].month,time_dt[0].day)
time_sec = np.array([(x-time_dt0).total_seconds() for x in time_dt])

time_dt = np.array(time_dt)

plt.figure(); 
plt.plot(time_sec/3600,var_1d_data['RemoveLongI2Cell'])
plt.grid(b=True)
plt.xlabel('time [h-UTC]')
plt.ylabel('Value')
plt.title('RemoveLongI2Cell')

master_time = np.arange(time_sec[0]-tres/2,time_sec[-1]+tres/2,tres)

#match_data_times(master_time,profiles,var_1d_data,time_sec)
prof_list = []
prof_names = []
if 'cross' in var_2d_data.keys():
    CrossPol = lp.LidarProfile(var_2d_data['cross'],time_sec,label='Cross Polarization Channel',descript = 'Cross Polarization\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL',StartDate=cal_start.date())
    CrossPol.time_resample(tedges=master_time,update=True,remainder=False)
    CrossPol.bg_subtract(BGIndex)
    prof_list = prof_list+[CrossPol]
    prof_names = prof_names+['cross']

if 'molecular' in var_2d_data.keys():
    Molecular = lp.LidarProfile(var_2d_data['molecular'],time_sec,label='Molecular Backscatter Channel',descript = 'Parallel Polarization\nMolecular Backscatter Returns',bin0=47,lidar='GV-HSRL',StartDate=cal_start.date())
    Molecular.time_resample(tedges=master_time,update=True,remainder=False)
    Molecular.bg_subtract(BGIndex)
    prof_list = prof_list+[Molecular]
    prof_names = prof_names+['molecular']

if 'combined_hi' in var_2d_data.keys():
    CombHi = lp.LidarProfile(var_2d_data['combined_hi'],time_sec,label='High Gain Total Backscatter Channel',descript = 'Parallel Polarization\nHigh Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL',StartDate=cal_start.date())
    CombHi.time_resample(tedges=master_time,update=True,remainder=False)
    CombHi.bg_subtract(BGIndex)
    prof_list = prof_list+[CombHi]
    prof_names = prof_names+['combined_hi']

if 'combined_lo' in var_2d_data.keys():       
    CombLo = lp.LidarProfile(var_2d_data['combined_lo'],time_sec,label='Low Gain Total Backscatter Channel',descript = 'Parallel Polarization\nLow Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL',StartDate=cal_start.date())            
    CombLo.time_resample(tedges=master_time,update=True,remainder=False)
    CombLo.bg_subtract(BGIndex)
    prof_list = prof_list+[CombLo]
    prof_names = prof_names+['combined_lo']

profs = dict(zip(prof_names,prof_list))

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
    
#fig_data = lp.pcolor_profiles([CombHi,Molecular],ylimits=[0,12],climits=[[1e-4,1e4],[1e-4,1e4]])  
#fig_data[1][1].plot(time_sec/3600,var_1d_data['RemoveLongI2Cell']/25,'b--')      
plt.show(block=False)


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
            i1 = np.concatenate((i1,np.array([-1])))
    
    print('%d.)  custom'%(ai+1))    
    cal_index = np.int(input('Select Interval (invalid number quits cal): '))
    
        
    
    if cal_index > len(i0):
        RunCal = False
    elif cal_index == len(i0):
        t_input1 = np.float(input('Enter start time in h-UTC: '))
        t_input2 = np.float(input('Enter stop time in h-UTC:'))
        time_range = [t_input1, t_input2]
        i0 = np.array(list(i0) + np.argmin(np.abs(time_sec-t_input1)))
        i1 = np.array(list(i1) + np.argmin(np.abs(time_sec-t_input2)))
    else:
        time_range = [5*60+time_sec[i0[cal_index]],time_sec[i1[cal_index]]-5*60]
if RunCal:
    for pname in profs.keys():
        profs[pname].slice_time(time_range)
        profs[pname].time_integrate()
    
    pHi = profs['combined_hi'].copy()   
    pLo = profs['combined_lo'].copy() 
    pHiLo = profs['combined_hi'].copy()
    
    pHi.divide_prof(profs['molecular'])   
    pLo.divide_prof(profs['molecular'])  
    pHiLo.divide_prof(profs['combined_lo'])  
  
    hi_smooth = savitzky_golay(pHi.profile.flatten(), 11, 5, deriv=0)  
    #    dhi_smooth = savitzky_golay(pHi.profile.flatten(), 11, 4, deriv=1)
        
    plt.figure()
    plt.plot(pHi.profile.flatten())
    plt.plot(hi_smooth)
    plt.title('Hi/Mol Ratio')
    plt.show(block=False)
    
    i_norm = 1000
    i_const = np.int(input('Make constant above index (e.g. 500): '))
    i_const_max = np.nonzero(pHi.SNR().flatten() < 0.2*pHi.SNR().flatten()[i_const])[0]
    i_const_max = i_const_max[np.nonzero(i_const_max > i_const)[0][0]]
    
    
    hi_diff_geo = hi_smooth
#    hi_diff_geo[i_const:] = hi_diff_geo[i_const]
#    plt.plot(hi_diff_geo)
    
    hi_diff_geo = fit_high_range(hi_diff_geo,pHi.profile_variance,i_const,i_const_max)
    plt.plot(hi_diff_geo,'--')
    
    hi_diff_geo = hi_diff_geo/hi_diff_geo[i_norm]
#    xfit = np.arange(hi_smooth[i_const:i_const_max].size)
#    yfit = hi_smooth[i_const:i_const_max]
#    wfit = 1.0/np.sqrt(pHi.profile_variance[0,i_const:i_const_max].flatten())
#    wfit[0:5] = 10*np.max(wfit)
#    #pfit = lp.polyfit_with_fixed_points(1,xfit,yfit, np.array([0]) ,np.array([geo_prof[200]]))
#    pfit = np.polyfit(xfit,yfit,1,w=wfit)
#    xprof = np.arange(hi_smooth[i_const:].size)
#    hi_diff_geo[i_const:] = np.polyval(pfit,xprof)
    
    lo_smooth = savitzky_golay(pLo.profile.flatten(), 11, 5, deriv=0) 
    
    plt.figure()
    plt.plot(pLo.profile.flatten())
    plt.plot(lo_smooth)
    plt.title('Lo/Mol Ratio')
    plt.show(block=False)
    
    i_const = np.int(input('Make constant above index (e.g. 100): '))
    i_const_max = np.nonzero(pLo.SNR().flatten() < 0.2*pLo.SNR().flatten()[i_const])[0]
    i_const_max = i_const_max[np.nonzero(i_const_max > i_const)[0][0]]
    
    
    lo_diff_geo = lo_smooth
#    lo_diff_geo[i_const:] = lo_diff_geo[i_const]
#    plt.plot(lo_diff_geo)
    
    lo_diff_geo = fit_high_range(lo_diff_geo,pLo.profile_variance,i_const,i_const_max)
    plt.plot(lo_diff_geo,'--')
    plt.show(block=False)
    
    lo_diff_geo = lo_diff_geo/lo_diff_geo[i_norm]
    
    save_cal = input("Save Calibrations[y/n]")
    
    if save_cal == 'y' or save_cal == 'Y':
        write_data = np.ones((hi_diff_geo.size,4))
        write_data[:,0] = np.arange(hi_diff_geo.size)
        write_data[:,1] = hi_diff_geo
        write_data[:,2] = lo_diff_geo
        
        header_str = 'profile data from ' + time_dt[i0[cal_index]].strftime('%d-%b-%y %H:%M')+' -->'+time_dt[i1[cal_index]].strftime('%H:%M UTC') +'\n'
        header_str = header_str+'Gains normalized range bin = %d'%i_norm + '\n'
        header_str = header_str+'bin #\tdiff_hi/mol\tdiff_lo/mol\tdiff_i2a/mol'
    
        save_cal_file = 'diff_geofile_'+time_dt[i0[cal_index]].strftime('%Y%m%dT%H%M')+'.geo'
#        save_file_path = '/Users/mhayman/Documents/Python/Lidar/'
        save_file_path = '/h/eol/mhayman/HSRL/hsrl_processing/hsrl_configuration/projDir/calfiles/'
        
        print('Saving to:\n  '+save_file_path+save_cal_file)        
        
        np.savetxt(save_file_path+save_cal_file,write_data,fmt='\t%d\t%.4e\t%.4e\t%.4e',header=header_str)
#    plt.figure()
#    plt.plot((profs['molecular'].profile/profs['combined_hi'].profile).flatten())
#    plt.title('molecular/combined hi')
#    
#    plt.figure(); 
#    plt.plot((profs['combined_lo'].profile/profs['combined_hi'].profile).flatten())
#    plt.title('combined lo/hi')
    
#    plt.show(block=False)
    
    
    
# check for power stability
# add 5 min buffers on either side
# trim profiles to the selected interval

# perform diff_geo

