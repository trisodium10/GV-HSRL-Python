# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 10:22:19 2017

@author: mhayman
"""

import numpy as np
import netCDF4 as nc4
import LidarProfileFunctions as lp
import datetime
import glob

import ptv.hsrl.denoise as denoise
from ptv.estimators.poissonnoise import poissonmodel0

import scipy.optimize

import matplotlib.pyplot as plt



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
    xprof = np.arange(profile[i_min:].size)
    profile_out = profile.copy()
    profile_out[i_min:] = np.polyval(pfit,xprof)    
    return profile_out
    
    
def load_raw_data(start_time,stop_time,var_2d_list,var_1d_list,basepath = '/scr/eldora1/rsfdata/hsrl/raw/',verbose=True,as_prof=True,loadQWP='fixed',date_reference=0):
    """
    loads GVHSRL raw data from netcdf files stored in basepath
    accepts a list of the variables to be loaded (their netcdf names)
    and returns a dict() of where the variable names correspond to the 
    data.
    
    start_time - datetime object describing the desired start time and date
    stop_time - datetime object describing the desired stop time and date
    
    If verbose is set to True the function reports the QWP status of each
    file.
    
    If as_prof = True, 2D data is returned as lidar profiles
    
    loadQWP - 'all':  load all data regardless of QWP status
              'fixed': load only fixed QWP data
              'rotating': load only rotating QWP data
    
    date_reference - datetime object reference for time axis.  this is needed when processing
        flights in chunks across midnight
    """
    
    var_1d_data = dict(zip(var_1d_list,[np.array([])]*len(var_1d_list)))
    var_2d_data = dict(zip(var_2d_list,[np.array([])]*len(var_2d_list)))
    timeD = np.array([])
    
    start_time_hr = datetime.datetime(year=start_time.year,month=start_time.month,day=start_time.day,hour=start_time.hour)
    if hasattr(date_reference,'year'):
        start_date = date_reference
    else:
        start_date = datetime.datetime(year=start_time.year,month=start_time.month,day=start_time.day,hour=0)   
    # build list of files falling on the requested dates
    day_iter =  start_time.date()
    SubFiles = []
    FileDate = []
    while day_iter <= stop_time.date():
        FilePath0 = basepath+day_iter.strftime('%Y/%m/%d/raw/')
        SubList0 = sorted(glob.glob(FilePath0+'*data*.nc'))
        SubFiles = SubFiles+SubList0                  # make a list of files falling on the dates corresponding to the search period
        FileDate = FileDate+len(SubList0)*[day_iter]  # collect dates associated with each file (it's just easier)
        day_iter = day_iter + datetime.timedelta(days=1)
    
    # iterate through files to determine if we should load them (based on time)
    # load them if we should
    for idir in range(len(SubFiles)):
        # get the file start time from the file name
        iTime = np.char.find(SubFiles[idir],'T',len(basepath),-1)
        Hour = np.double(SubFiles[idir][iTime+1:iTime+3])
        File_Time = datetime.datetime(*FileDate[idir].timetuple()[:4])+datetime.timedelta(hours=Hour)
        # check if this file is in the search time period
        if File_Time >= start_time_hr and File_Time <= stop_time:
            if verbose:
                print( 'Loading %d UT' %Hour)
                print( '['+SubFiles[idir]+']')
                
            f = nc4.Dataset(SubFiles[idir],'r')
            
            # System for time array still needs work
            timeD0 = np.array(f.variables['DATA_time'][:]).astype(np.int)
            timeD0[:,6] = timeD0[:,6]*1e3+timeD0[:,7]
            time_dt0 = [datetime.datetime(*x) for x in timeD0[:,:7]]
            time_sec0 = np.array([(x-time_dt0[0]).total_seconds() for x in time_dt0])
            
            if any('polarization' in s for s in f.variables):
                QWP = f.variables['polarization'][:].copy()               # QWP rotation angle
                # Calculate the QWP rotation rate to determine its status
                meanQWPrate = np.median(np.diff(QWP)/np.diff(time_sec0)*180/np.pi)
                meanQWPvalue = np.median(QWP*180/np.pi)
                if meanQWPrate >= 10:
                    QWP_Status = 'rotating'
                else:
                    QWP_Status = 'fixed'
                if verbose:
                    print( 'QWP Rotation Rate: %f deg/sec' %meanQWPrate)
                    print( 'Average QWP Position: %f deg' %meanQWPvalue)
            else:
                QWP_Status = 'fixed'
                
            if loadQWP == 'fixed' and QWP_Status == 'rotating':
                load_set = False
            elif loadQWP == 'rotating' and QWP_Status == 'fixed':
                load_set = False
            else:
                load_set = True
            
            
            if load_set:
            
                for var in var_1d_data.keys():
                    if any(var in s for s in f.variables):
                        var_data = lp.ncvar(f,var)
                        if len(var_1d_data[var]) > 0:
                            var_1d_data[var] = np.concatenate((var_1d_data[var],var_data)) 
                        else:
                            var_1d_data[var] = var_data.copy()
                    
                
                
                if len(timeD)>0:
                    timeD = np.vstack((timeD,np.array(f.variables['DATA_time'][:].astype(np.int))));          # time array [year,month,day,hour,minute,second,msec,usec]
                else:
                    timeD = np.array(f.variables['DATA_time'][:]).astype(np.int)
    
                
                
                if verbose:
                    print( 'Profile Integration Time: %f seconds' %np.median(np.diff(time_sec0)))
                    print( 'Processing QWP as %s\n' %QWP_Status)
                    
                
                for var in var_2d_data.keys():
                    if len(var_2d_data[var]) == 0:
                        var_2d_data[var] = f.variables[var][:].copy()
                    else:
                        var_2d_data[var] = np.vstack((var_2d_data[var],f.variables[var][:]))
            elif verbose:
                print( 'Skipping %d UT' %Hour)
                print( '['+SubFiles[idir]+']')
                print( 'Profile Integration Time: %f seconds' %np.median(np.diff(time_sec0)))
                print( 'Processing QWP as %s\n' %QWP_Status)
                    
        
            f.close()
    if len(timeD) > 0:
        #timeD = timeD[:,:8]
        timeD[:,6] = timeD[:,6]*1e3+timeD[:,7]
        time_dt = [datetime.datetime(*x) for x in timeD[:,:7]]
        time_dt0 = start_date #datetime.datetime(time_dt[0].year,time_dt[0].month,time_dt[0].day)
        time_sec = np.array([(x-time_dt0).total_seconds() for x in time_dt])
        
        time_dt = np.array(time_dt)
        
        # trim the data to match the requested times
        ifirst = np.nonzero(time_dt < start_time)[0]
        if len(ifirst) == 0:
            ifirst = 0
        else:
            ifirst = ifirst[-1]
        
        
        ilast = np.nonzero(time_dt > stop_time)[0]
        if len(ilast) == 0:
            ilast = -1
        else:
            ilast = ilast[0]
            
        # trim 2d data
        for var in var_2d_data.keys():
            var_2d_data[var] = var_2d_data[var][ifirst:ilast,:]
            
        # trim 1d data
        for var in var_1d_data.keys():
            var_1d_data[var] = var_1d_data[var][ifirst:ilast]
            
        # trim time arrays
        timeD = timeD[ifirst:ilast,:]
        time_dt = time_dt[ifirst:ilast]
        time_sec = time_sec[ifirst:ilast]
        
        
        time_list = [timeD,time_dt,time_sec]
    
    
    
        # if requested, load 2D data as data types LidarProfile
        if as_prof:
            prof_list = []
            prof_names = []
            
            for var in var_2d_data.keys():
                # check for known channel names to fill in the lidar profile definitions
                if var == 'molecular':
                    plabel = 'Molecular Backscatter Channel'
                    pdescript = 'Parallel Polarization\nMolecular Backscatter Returns'
                elif var == 'cross':
                    plabel = 'Cross Polarization Channel'
                    pdescript = 'Cross Polarization\nCombined Aerosol and Molecular Returns'
                elif var == 'combined_hi':
                    plabel = 'High Gain Total Backscatter Channel'
                    pdescript = 'Parallel Polarization\nHigh Gain\nCombined Aerosol and Molecular Returns'
                elif var == 'combined_lo':
                    plabel = 'Low Gain Total Backscatter Channel'
                    pdescript = 'Parallel Polarization\nLow Gain\nCombined Aerosol and Molecular Returns'
                else:
                    plabel = 'Unassigned Channel'
                    pdescript = 'Unassigned Channel Description'
                
                prof_list = prof_list+[lp.LidarProfile(var_2d_data[var],\
                    time_sec,label=plabel,\
                    descript = pdescript,\
                    bin0=37,lidar='GV-HSRL',StartDate=start_time.date())]
                prof_names = prof_names + [var]
                
            profs_2d = dict(zip(prof_names,prof_list))
        else:
            profs_2d = var_2d_data
    else:
        time_list = []
        print('No lidar data found in')
        print('   ' + basepath)
        print('   ' + 'last search dir: ' + FilePath0)
        
        # return empties if no data
        time_list = []
        var_1d_data = {}
        profs_2d = {}
    
    return time_list, var_1d_data, profs_2d

def load_aircraft_data(filename,var_list):
    """
    Loads data from aircraft netcdf
    filename - string with path and filename to load
    var_list - strings of the netcdf variables to load
    """    

    print('Loading aircraft data file:')
    print('   '+filename)
    
    var_data = dict(zip(var_list,[np.array([])]*len(var_list)))
    
    try:
        f = nc4.Dataset(filename,'r')
        
        i_nan = []  # list of invalid values from the data system     
        for var in var_data.keys():
            if any(var in s for s in f.variables):
                data = lp.ncvar(f,var)
                if len(var_data[var]) > 0:
                    var_data[var] = np.concatenate((var_data[var],data)) 
                else:
                    var_data[var] = data.copy()
                i_nan.extend(list(np.nonzero(var_data[var]==-32767)[0]))
        i_nan = np.unique(np.array(i_nan))  # force list to be unique
        
        # delete data where nans are present
        for var in var_data.keys():
            var_data[var] = np.delete(var_data[var],i_nan)
        print('Aircraft Time Data: %s' %f.variables['Time'].units)
        time_ref = datetime.datetime.strptime(f.variables['Time'].units[:],'seconds since %Y-%m-%d %H:%M:%S +0000')
        f.close()
    except RuntimeError:
        print('Aircraft data file NOT found')
        time_ref=[]
    
    return var_data,time_ref

def interp_aircraft_data(master_time,aircraft_data):
    """
    Interpolates aircraft data to the lidar time grid
    master_time - the lidar time grid in seconds since the processing start date at 0:00 UTC
    aircraft_data - dict of variables loaded from the aircraft data file including
        'Time' variable
    """
    air_data_new = {}
    for var in aircraft_data.keys():
        if var != 'Time':
            air_data_new[var] = np.interp(master_time,aircraft_data['Time'],aircraft_data[var])
    air_data_new['Time'] = master_time.copy()  # update the time variable
    return air_data_new

def var_time_resample(tedges,var_time,varlist,average=True,remainder=False):
    """
    regrids 1D data onto a predefined time grid
    tedges - desired time grid
    var_time - time array of current variables
    varlist - list or dict of the variables to be regridded
    average - (if True) bin by averaging, (if False) bin by summing
    remainder - (if True) returns remainder that didn't fit onto the grid
    """    
    
    varlist_new = []
    remlist = []
    var_names = []
    # deterimine the bins for each time index       
    itime = np.digitize(var_time,tedges)
    # Only run if the profiles fit in the master timing array (tedges), otherwise everything is a remainder
    if np.sum((itime>0)*(itime<tedges.size))!=0:
        
#                iremain = np.nonzero(self.time > tedges[-1])[0]
        
#        iremain = np.int(np.max(itime))  # the remainder starts at the maximum bin where data exists
        iremain = tedges.size
#                if not remainder and iremain < tedges.size:
#                    iremain = iremain+1  # make sure to include the last bin of data if remainder isn't being used.
        iremainList = np.nonzero(var_time > tedges[iremain-1])[0]
        
#        iprofstart = np.int(np.max(np.array([1,np.min(itime)])))
        iprofstart = 1
#                print('start index: %d\nstop index: %d'%(iprofstart,iremain))
        
#                profNew = np.zeros((np.size(tedges)-1,self.profile.shape[1]))
#                timeNew = 0.5*tedges[1:]+0.5*tedges[:-1] 
        timeNew = -np.diff(tedges[iprofstart-1:iremain])*0.5+tedges[iprofstart:iremain]
        for xi,var in enumerate(varlist):
            if isinstance(varlist,dict):
                varOld = varlist[var].copy()
            else:
                varOld = var.copy()
            
            if varOld.ndim == 1:
                varOld = varOld[:,np.newaxis]
            
            varNew = np.zeros((iremain-iprofstart,varOld.shape[1]))
            
##                itimeNew = np.arange(iprofstart,iremain)
#            var_profNew = np.zeros(profNew.shape)
#            shot_countNew = np.zeros(timeNew.shape)
#            self.NumProfList = np.zeros(timeNew.shape)
                  

            for ai in range(np.size(timeNew)):
                if hasattr(varOld,'mask'):
    #                        NumProf = np.nansum(np.logical_not(self.profile[itime == ai+iprofstart,:].mask),axis=0)
    #                        NumProfDiv = NumProf.copy()
    #                        NumProfDiv[np.nonzero(NumProf==0)] = 1
                    NumProf = np.nanmax(np.nansum(np.logical_not(varOld[itime == ai+iprofstart,:].mask),axis=0))
                    
                    if NumProf == 0:
                        NumProfDiv = 1.0
                    else:
                        NumProfDiv = NumProf
                else:
                    NumProf = varOld[itime == ai+iprofstart,:].shape[0]
                    if NumProf == 0:
                        NumProfDiv = 1.0
                    else:
                        NumProfDiv = np.float(NumProf)
                if not average:
                    NumProfDiv = 1.0
    
                varNew[ai,:] = np.nansum(varOld[itime == ai+iprofstart,:],axis=0)/NumProfDiv
                
            if remainder:
                varRem = varOld[iremainList,:].copy()
                if varRem.shape[1] == 1:
                    remlist = remlist + [varRem.flatten()]
                else:
                    remlist = remlist + [varRem]
            
            if varNew.shape[1] == 1:
                varlist_new = varlist_new + [varNew.flatten()]
            else:
                varlist_new = varlist_new + [varNew]
            
            if isinstance(varlist,dict):
                var_names = var_names + [var]
        
        if isinstance(varlist,dict):
            var_return = dict(zip(var_names,varlist_new))  
        else:
            var_return = varlist_new
                
        if remainder:
            return timeNew,var_return,varRem
        else:
            return timeNew,var_return

#def var_time_resample(tedges,var_time,varlist,average=True,remainder=False):
#    """
#    regrids 1D data onto a predefined time grid
#    tedges - desired time grid
#    var_time - time array of current variables
#    varlist - list or dict of the variables to be regridded
#    average - (if True) bin by averaging, (if False) bin by summing
#    remainder - (if True) returns remainder that didn't fit onto the grid
#    """    
#    
#    varlist_new = []
#    remlist = []
#    var_names = []
#    # deterimine the bins for each time index       
#    itime = np.digitize(var_time,tedges)
#    # Only run if the profiles fit in the master timing array (tedges), otherwise everything is a remainder
#    if np.sum((itime>0)*(itime<tedges.size))!=0:
#        
##                iremain = np.nonzero(self.time > tedges[-1])[0]
#        iremain = np.int(np.max(itime))  # the remainder starts at the maximum bin where data exists
##                if not remainder and iremain < tedges.size:
##                    iremain = iremain+1  # make sure to include the last bin of data if remainder isn't being used.
#        iremainList = np.nonzero(var_time > tedges[iremain-1])[0]
#        iprofstart = np.int(np.max(np.array([1,np.min(itime)])))
##                print('start index: %d\nstop index: %d'%(iprofstart,iremain))
#        
##                profNew = np.zeros((np.size(tedges)-1,self.profile.shape[1]))
##                timeNew = 0.5*tedges[1:]+0.5*tedges[:-1] 
#        timeNew = -np.diff(tedges[iprofstart-1:iremain])*0.5+tedges[iprofstart:iremain]
#        for xi,var in enumerate(varlist):
#            if isinstance(varlist,dict):
#                varOld = varlist[var].copy()
#            else:
#                varOld = var.copy()
#            
#            if varOld.ndim == 1:
#                varOld = varOld[:,np.newaxis]
#            
#            varNew = np.zeros((iremain-iprofstart,varOld.shape[1]))
#            
###                itimeNew = np.arange(iprofstart,iremain)
##            var_profNew = np.zeros(profNew.shape)
##            shot_countNew = np.zeros(timeNew.shape)
##            self.NumProfList = np.zeros(timeNew.shape)
#                  
#
#            for ai in range(np.size(timeNew)):
#                if hasattr(varOld,'mask'):
#    #                        NumProf = np.nansum(np.logical_not(self.profile[itime == ai+iprofstart,:].mask),axis=0)
#    #                        NumProfDiv = NumProf.copy()
#    #                        NumProfDiv[np.nonzero(NumProf==0)] = 1
#                    NumProf = np.nanmax(np.nansum(np.logical_not(varOld[itime == ai+iprofstart,:].mask),axis=0))
#                    
#                    if NumProf == 0:
#                        NumProfDiv = 1.0
#                    else:
#                        NumProfDiv = NumProf
#                else:
#                    NumProf = varOld[itime == ai+iprofstart,:].shape[0]
#                    if NumProf == 0:
#                        NumProfDiv = 1.0
#                    else:
#                        NumProfDiv = np.float(NumProf)
#                if not average:
#                    NumProfDiv = 1.0
#    
#                varNew[ai,:] = np.nansum(varOld[itime == ai+iprofstart,:],axis=0)/NumProfDiv
#                
#            if remainder:
#                varRem = varOld[iremainList,:].copy()
#                if varRem.shape[1] == 1:
#                    remlist = remlist + [varRem.flatten()]
#                else:
#                    remlist = remlist + [varRem]
#            
#            if varNew.shape[1] == 1:
#                varlist_new = varlist_new + [varNew.flatten()]
#            else:
#                varlist_new = varlist_new + [varNew]
#            
#            if isinstance(varlist,dict):
#                var_names = var_names + [var]
#        
#        if isinstance(varlist,dict):
#            var_return = dict(zip(var_names,varlist_new))  
#        else:
#            var_return = varlist_new
#                
#        if remainder:
#            return timeNew,var_return,varRem
#        else:
#            return timeNew,var_return
            
def get_TP_from_aircraft(air_data,profile,telescope_direction=[],lidar_tilt=[0,4.0]):
    """
    Returns a temperature and pressure profile corresponding to
    dimensions of profile
    The profile should already be converted to altitude coordinates from range
    
    The aircraft temperature and pressure measurements are contained
    in air_data
    
    If telescope_direction is passed in, assume this is range centered 
    processing
    
    lidar_tilt defines the lidar tilt angles relative to the aircraft.  It
    is only used in range centered processing.  It is defined 
    [forward/aft,port/starboard] in degrees
    """
    
    # make sure the time axes of the aircraft data aline with the 
    # lidar profile data
    
    if len(telescope_direction) == 0:
        # process in altitude centered configuration
        air_data_int = interp_aircraft_data(profile.time,air_data)
        aircraft_temp = air_data_int['ATX']+273.15  # convert aircraft temperature to K 
        b_T = (aircraft_temp + air_data_int['GGALT']*0.0065)[:,np.newaxis]
        TempAir = b_T-0.0065*profile.range_array[np.newaxis,:]
        PresAir = air_data_int['PSXC'][:,np.newaxis]*100*(aircraft_temp[:,np.newaxis]/TempAir)**(-5.5)  # pressure in Pa from PSXC in hPa
    else:
        # process in range centered configuration
        air_data_int = interp_aircraft_data(profile.time,air_data)
        aircraft_temp = air_data_int['ATX']+273.15  # convert aircraft temperature to K 
        b_T = (aircraft_temp + air_data_int['GGALT']*0.0065)[:,np.newaxis]        
        telescope_direction = np.sign(telescope_direction-0.5)

        # create a 2D array of all the raw altitude data caputured by the lidar
#        alt_raw = \
#            (profile.range_array[:,np.newaxis]*telescope_direction[np.newaxis,:]\
#            *np.cos((air_data_int['ROLL'][np.newaxis,:]+lidar_tilt[1])*np.pi/180)*np.cos((air_data_int['PITCH'][np.newaxis,:]+lidar_tilt[0])*np.pi/180)\
#            +air_data_int['GGALT'][np.newaxis,:]).T
        alt_raw = \
            (profile.range_array[:,np.newaxis]*telescope_direction[np.newaxis,:]\
            *np.cos((air_data_int['ROLL']-lidar_tilt[1]*telescope_direction)*np.pi/180)*np.cos((air_data_int['PITCH']+lidar_tilt[0])*np.pi/180)\
            +air_data_int['GGALT'][np.newaxis,:]).T

        TempAir = b_T-0.0065*alt_raw
        PresAir = air_data_int['PSXC'][:,np.newaxis]*100*(aircraft_temp[:,np.newaxis]/TempAir)**(-5.5)  # pressure in Pa from PSXC in hPa
    
    
    temp = profile.copy()
    temp.profile = TempAir.copy()
    temp.profile_variance = (temp.profile*0.1)**2  # force SNR of 10 in sonde profile.
    temp.ProcessingStatus = []     # status of highest level of lidar profile - updates at each processing stage
    temp.lidar = 'Aircraft measurement'
    
    temp.diff_geo_Refs = ['none']           # list containing the differential geo overlap reference sources (answers: differential to what?)
    temp.profile_type =  '$K$'
    
    temp.bg = np.zeros(temp.bg.shape) # profile background levels
    
    temp.descript = 'Ideal Atmosphere Temperature in K'
    temp.label = 'Temperature'
    
    pres = profile.copy()
    pres.profile = PresAir
    pres.profile_variance = (pres.profile*0.1)**2  # force SNR of 10 in sonde profile.
    pres.ProcessingStatus = []     # status of highest level of lidar profile - updates at each processing stage
    pres.lidar = 'Aircraft measurement'
    
    pres.diff_geo_Refs = ['none']           # list containing the differential geo overlap reference sources (answers: differential to what?)
    pres.profile_type =  '$Pa$'
    
    pres.bg = np.zeros(temp.bg.shape) # profile background levels
    
    pres.descript = 'Ideal Atmosphere Pressure in Pa'
    pres.label = 'Pressure'
    
    return temp, pres
    
def delete_indices(in_dict,indices):
    """
    deletes the indices of all elements of the dictionary
    """
    out_dict = {}
    for var in in_dict.keys():
        out_dict[var] = np.delete(in_dict[var],indices)
    return out_dict
    
    
def DenoiseMolecular(MolRaw,beta_m_sonde=np.array([np.nan]),
                     geo_data=dict(geo_prof=np.array([1])),
                    MaxAlt=np.nan,n=1,start_time=0,end_time=np.nan,
                    verbose=False,accel = False,tv_lim =[0.4, 1.8],N_tv_pts=48,
                    bg_index = -50,geo_key='geo_prof',MolGain_Adj = 0.75,
                    plot_result=False,eps_opt = 1e-5,tv_bound=True):
    """
    Use Willem Marais' functions to denoise the molecular signal in an 
    HSRL signal.
    MolRaw - raw profile
    beta_m_sonde - estimated molecular backscatter coefficient.  Not used if
                not provided.
    geo_data - geometric overlap function.  Not used if not provided
    MaxAlt - maximum altitude in meters to fit.  Defaults to full profile.
    n = number of time profiles to process at a time
    start_time - the profile time to start on.  Defaults to first profile.
    end_time - the profile time to end on.  Defaults to last profile
    verbose - set to True to have text output from optimizer routine
    accel - attempt to accelerate by using previous TV results
    tv_lim - limits for the tv search space.  default of [0.4, 1.8] are obtained
        for the DLB-HSRL  other settings may be desirable for other lidar
    N_tv_pts - number of tv points to evalute in the tv_lim space.  Defaults
        to 48 used in DLB-HSRL
    bg_index - index above which to treat as background
    geo_key - key to access geo overlap data
    MolGain_Adj - factor to multiply the molecular gain to get a decent agreement
        between estimated backscatter and observed.
    plot_result - plot results of each iteration (for debugging only)
    eps_opt - optimization precision.
    tv_bound - bound optimization limits to those given in tv_lim.  Only used
        if accel = True.  Otherwise the accelerated tv_lims are unbounded.
    """
    
#    bg_index = -50    
    if not np.isnan(end_time):
        MolRaw.slice_time([start_time,end_time])
    if any('Background Subtracted' in s for s in MolRaw.ProcessingStatus):
        BG_Sub_Flag = True
        MolRaw.profile = MolRaw.profile+MolRaw.bg[:,np.newaxis]
        Mol_BG_Total = MolRaw.bg*MolRaw.NumProfList
    else:
        BG_Sub_Flag = False    
        Mol_BG_Total = np.mean(MolRaw.profile[:,bg_index:],axis=1)*MolRaw.NumProfList  # get background before we remove the top of the profile
    
    
    
    
    if not np.isnan(MaxAlt):    
        MolRaw.slice_range(range_lim=[0,MaxAlt]) 

    
    
    
    MolDenoise = MolRaw.copy()
    MolDenoise.label = 'Denoised Molecular Backscatter Channel'
    MolDenoise.descript = 'Total Variation Denoised\nUnpolarization\nMolecular Backscatter Returns'
#    MolDenoise.profile_type = '$m^{-3}$'
        
    
    if n  > MolRaw.time.size:
        n = MolRaw.time.size
        
    tune_list = []
    
    fit_range_array = MolRaw.range_array.copy()
    fit_range_array[np.nonzero(fit_range_array==0)] = 1
    
    percent_index = 0
    percent_array = np.arange(0,110,5)
    
    for i_prof in range(np.ceil(MolRaw.time.size*1.0/n).astype(np.int)):
        
        percent_complete = 100.0*i_prof/np.ceil(MolRaw.time.size*1.0/n)
        if percent_complete >= percent_array[percent_index]:
            percent_index = percent_index+1
            print('%d %% complete'%percent_complete)
        
    
        istart = i_prof*n
        iend = np.min(np.array([istart + n,MolRaw.time.size]))

#        MolGain_Adj = 0.75 # 0.5

        MolFit = (MolRaw.profile[istart:iend,:]*MolRaw.NumProfList[istart:iend,np.newaxis]).astype (np.int)
        NumProf = MolRaw.NumProfList[istart:iend,np.newaxis]
        
        Mol_BG = Mol_BG_Total[istart:iend]# np.mean(MolFit[:,bg_index:],axis=1)


        # Check what a priori data the user is providing
        try:
            geo_len = len(geo_data[geo_key])
            if geo_len == 1:
                # no geo data provided.
                geo_est = np.ones(MolFit.shape[1])*geo_data[geo_key]
            else:
                if not 'Nprof' in geo_data.keys():
                    geo_data['Nprof'] = 1
                if geo_data[geo_key].ndim == 1:
                    # 1d array for geo overlap estimate
                    geofun = 1/np.interp(MolRaw.range_array,geo_data[geo_key][:,0],geo_data[geo_key][:,1])
                    geo_est = 1.0/geo_data['Nprof']*geofun
                else:
                    # 2d array for geo overlap estimate
                    # typically used for up and down (airborne) pointing
                    geofun0 = 1.0/geo_data[geo_key][istart:iend,:]
                    geofun = np.zeros(MolFit.shape)
                    for iprof in range(geofun.shape[0]):
                        geofun[iprof,:] = np.interp(MolRaw.range_array,geo_data['range_array'],geofun0[iprof,:])
                        
                    geo_est = 1.0/geo_data['Nprof']*geofun
    #            geo_est = MolRaw.mean_dt/geo_data['Nprof']/geo_data['tres']*geofun
        except TypeError:
            geo_est = np.ones(MolFit.shape[1])*geo_data
            
        if hasattr(beta_m_sonde,'profile'):
            # Estimated molecular backscatter is passed to the function
#            if not np.isnan(MaxAlt):    
#                beta_m_sonde.slice_range(range_lim=[0,MaxAlt]) 
            beta_m_2D0 = beta_m_sonde.profile[istart:iend,:]
            beta_m_2D = np.zeros(MolFit.shape)
            for iprof in range(beta_m_2D0.shape[0]):
                beta_m_2D[iprof,:] = np.interp(MolRaw.range_array,beta_m_sonde.range_array,beta_m_2D0[iprof,:])
            
        else:
            beta_m_2D = np.ones(MolFit.shape)
        
#        print(beta_m_2D.shape)
#        print(geo_est.shape)        
#        print(fit_range_array.shape)
        
        
        # Create the Poisson thin object so that we can do cross-validation
        poisson_thn_obj = denoise.poissonthin (MolFit.T, p_trn_flt = 0.5, p_vld_flt = 0.5)
        
        # define coefficients in fit
        A_arr = (NumProf*MolGain_Adj*beta_m_2D*geo_est/fit_range_array**2).T
        A_arr[np.nonzero(A_arr==0)] = 1e-20
        A_arr[np.nonzero(np.isnan(A_arr))] = 1e-20
        A_arr[np.nonzero(np.isinf(A_arr))] = 1e-20
        
#        # for debugging
        if plot_result:
            import matplotlib.pyplot as plt
            plt.figure()
            plt.semilogy(MolFit.flatten())
            plt.semilogy(A_arr.flatten()+Mol_BG[0])
        

        sparsa_cfg_obj = denoise.sparsaconf (eps_flt = eps_opt, verbose_int = 1e6) # change eps from 1e-5 to 1e-7
        
        # check if the fit data is 1D or 2D.  1D can be run faster.
        if MolFit.shape[0] == 1:
            # Use for 1D denoising
            est_obj = poissonmodel0 (poisson_thn_obj, A_arr = A_arr, b_arr=Mol_BG[np.newaxis], log_model_bl = True, penalty_str = 'condatTV', 
                sparsaconf_obj = sparsa_cfg_obj)
        else:
            # Use for 2D denoising
            est_obj = poissonmodel0 (poisson_thn_obj, A_arr = A_arr, b_arr=Mol_BG[np.newaxis], log_model_bl = True, penalty_str = 'TV', 
                sparsaconf_obj = sparsa_cfg_obj)    
            
        # Create the denoiser object
        if MolFit.shape[0] == 1:
            # Use for 1D denoising:  
            # Defaults: log10_reg_lst = [-2, 2], nr_reg_int = 48
        
            if accel and i_prof > 0:
                if tv_bound:
                    tv_min = np.max([log_tune[np.argmin(valid_val)]*0.8,tv_lim[0]])
                    tv_max = np.min([log_tune[np.argmin(valid_val)]*1.2,tv_lim[1]])
                    tv_reg = [tv_min, log_tune[np.argmin(valid_val)]*1.2]
                else:
                    tv_reg = [log_tune[np.argmin(valid_val)]*0.8, tv_max]  # log of range of TV values to test
                nr_int = 5  # number of TV values to test
            else:
                tv_reg = tv_lim  # log of range of TV values to test
#                nr_int = 17  # number of TV values to test
                if tv_lim[0] == tv_lim[1]:
                    nr_int = 1
                else:
                    nr_int = N_tv_pts
                
            denoise_cnf_obj = denoise.denoiseconf (log10_reg_lst = tv_reg, nr_reg_int =nr_int, 
                pen_type_str = 'condatTV', verbose_bl = verbose)
        else:
            # Use for 2D denoising
            denoise_cnf_obj = denoise.denoiseconf (log10_reg_lst = [-2, 2], nr_reg_int = 48, 
                pen_type_str = 'TV', verbose_bl = verbose)
            
        denoiser_obj = denoise.denoisepoisson (est_obj, denoise_cnf_obj)
        # Start the denoising
        denoiser_obj.denoise ()
        
        MolDenoise.profile[istart:iend,:] = denoiser_obj.getdenoised().T/NumProf      
        
        
        log_tune,valid_val = denoiser_obj.get_validation_loss()
        tune_list.extend([[log_tune,valid_val]])  # store the results from the tuning parameters
        if plot_result:        
            plt.figure()
            plt.semilogy(MolFit.flatten())
            plt.semilogy(denoiser_obj.getdenoised().flatten(),'--')
            
            plt.figure()
            plt.plot(log_tune,valid_val)
            plt.show()
    
    if BG_Sub_Flag:
        MolRaw.profile = MolRaw.profile - MolRaw.bg[:,np.newaxis]
    
    MolDenoise.bg = Mol_BG_Total/MolRaw.NumProfList
    MolDenoise.bg_var = Mol_BG_Total/MolRaw.NumProfList**2
    MolDenoise.profile = MolDenoise.profile-MolDenoise.bg[:,np.newaxis]
    
    MolDenoise.ProcessingStatus.extend(['Applied range PTV denoising'])
    
    return MolDenoise,tune_list




def DenoiseTime(ProfRaw,MinAlt=0,MaxAlt=1e8,n=1,start_time=-1,end_time=1e16,
                    verbose=False,accel = False,tv_lim =[0.4, 3.0],N_tv_pts=48,
                    eps_opt = 1e-5,plot_result=False):
    """
    Use Willem Marais' functions to denoise the molecular signal in an 
    HSRL signal.
    MolRaw - raw profile
    beta_m_sonde - estimated molecular backscatter coefficient.  Not used if
                not provided.
    geo_data - geometric overlap function.  Not used if not provided
    MaxAlt - maximum altitude in meters to fit.  Defaults to full profile.
    MinAlt - minimum altitude in meters to fit.  Defaults to first gate.
    n = number of time profiles to process at a time
    start_time - the profile time to start on.  Defaults to first profile.
    end_time - the profile time to end on.  Defaults to last profile
    verbose - set to True to have text output from optimizer routine
    accel - attempt to accelerate by using previous TV results
    tv_lim - limits for the tv search space.  default of [0.4, 3.0] are obtained
        for the WV-DIAL  other settings may be desirable for other lidar
    N_tv_pts - number of tv points to evalute in the tv_lim space.  Defaults
        to 48 used in WV-DIAL
    """ 

    ProfDenoise = ProfRaw.copy()
#    ProfDenoise.slice_time([start_time,end_time])
#    ProfDenoise.slice_range(range_lim=[0,MaxAlt]) 
    range_index_end = np.argmin(np.abs(MaxAlt-ProfDenoise.range_array))
#    range_index_end = np.ceil(range_index_end-range_index_start*1.0/n).astype(np.int)
    range_index_start = np.argmin(np.abs(MinAlt-ProfDenoise.range_array))
    ProfDenoise.label = 'Denoised ' + ProfRaw.label
    ProfDenoise.descript = 'Total Variation Denoised in Time\n' + ProfRaw.descript
#    MolDenoise.profile_type = '$m^{-3}$'
        
    # check if background has been subtracted.  If so we need to add it back in
    # before denoising.
    if any('Background Subtracted' in s for s in ProfRaw.ProcessingStatus):
        BG_Sub_Flag = True
        ProfDenoise.profile = ProfDenoise.profile+ProfDenoise.bg[:,np.newaxis]
        ProfBG = ProfDenoise.bg[:,np.newaxis]
    else:
        BG_Sub_Flag = False
        ProfBG = np.mean(ProfDenoise.profile[:,-200:],axis=1)[:,np.newaxis]
    
    if n  > ProfRaw.range_array.size:
        n = ProfRaw.range_array.size
        
    tune_list = []
    
    report_percents = np.arange(0,105,5)
    report_index = 0    
    
    for i_prof in range(range_index_start,range_index_end+n,n):
#    for i_prof in range(np.ceil(ProfDenoise.range_array.size*1.0/n).astype(np.int)):
    
#        istart = i_prof*n
#        iend = np.min(np.array([istart + n,ProfRaw.range_array.size]))
        istart = i_prof
        iend = np.min(np.array([istart + n,ProfRaw.range_array.size]))
        
        percent_complete = (i_prof-range_index_start)*100.0/(range_index_end-range_index_start)
#        print('%f%% complete'%percent_complete)
        if percent_complete >= report_percents[report_index]:
            print('%f%% complete'%percent_complete)
            report_index = report_index + 1
        

#        if BG_Sub_Flag:
#            ProfFit = ((ProfRaw.profile[:,istart:iend]+ProfRaw.bg[:,np.newaxis])*ProfRaw.NumProfList[:,np.newaxis]).astype (np.int)
#        else:
#            ProfFit = (ProfRaw.profile[:,istart:iend]*ProfRaw.NumProfList[:,np.newaxis]).astype (np.int)
        if hasattr(ProfDenoise.profile,'mask'):
            ProfFit = (ProfDenoise.profile.data[:,istart:iend]*ProfDenoise.NumProfList[:,np.newaxis]).astype (np.int)
        else:
            ProfFit = (ProfDenoise.profile[:,istart:iend]*ProfDenoise.NumProfList[:,np.newaxis]).astype (np.int)
        NumProf = ProfDenoise.NumProfList[:,np.newaxis]
#        NumProf = ProfDenoise.NumProfList[istart:iend,np.newaxis]

        ProfFit[np.nonzero(ProfFit < 0)] = 0  # don't allow negative numbers

        
        # Create the Poisson thin object so that we can do cross-validation
        poisson_thn_obj = denoise.poissonthin (ProfFit, p_trn_flt = 0.5, p_vld_flt = 0.5)
        
        # define coefficients in fit
        A_arr = NumProf.astype(np.double) #ProfFit.astype(np.double) #NumProf
        A_arr[np.nonzero(A_arr==0)] = 1e-20
        A_arr[np.nonzero(np.isnan(A_arr))] = 1e-20
        A_arr[np.nonzero(np.isinf(A_arr))] = 1e-20

        sparsa_cfg_obj = denoise.sparsaconf (eps_flt = eps_opt, verbose_int = 1e6)
        
        # check if the fit data is 1D or 2D.  1D can be run faster.
        if ProfFit.shape[1] == 1:
            # Use for 1D denoising
            est_obj = poissonmodel0 (poisson_thn_obj, A_arr = A_arr, b_arr=ProfBG, log_model_bl = True, penalty_str = 'condatTV', 
                sparsaconf_obj = sparsa_cfg_obj)    # b_arr=Mol_BG[np.newaxis]
        else:
            # Use for 2D denoising
            est_obj = poissonmodel0 (poisson_thn_obj, A_arr = A_arr, log_model_bl = True, penalty_str = 'TV', 
                sparsaconf_obj = sparsa_cfg_obj)    #, b_arr=Mol_BG[np.newaxis]
                
        # while loop for accelerated denoising -
        # if end points are the selected values, adjust the limits and try again
        try_denoise = True
        accel_iter = 0  # number of times through the accelerated denoising loop
        max_iter = 10  # maximum number of iterations allowed through the denoising loop
        while try_denoise:        
            
            # Create the denoiser object
            if ProfFit.shape[1] == 1:
                # Use for 1D denoising:  
                # Defaults: log10_reg_lst = [-2, 2], nr_reg_int = 48
            
                if accel and i_prof > range_index_start:      
                    if verbose:
                        print('iteration: %d'%accel_iter)
                    
                    tv_reg = [log_tune[np.argmin(valid_val)]*0.8, log_tune[np.argmin(valid_val)]*1.2]  # log of range of TV values to test
                    nr_int = 5  # number of TV values to test
                    adapt_adj = True
                    accel_iter = accel_iter + 1                    
                else:
                    tv_reg = tv_lim  # log of range of TV values to test
                    nr_int = N_tv_pts   # number of TV values to test
                    adapt_adj = False
                    try_denoise = False
                    
                denoise_cnf_obj = denoise.denoiseconf (log10_reg_lst = tv_reg, nr_reg_int =nr_int, 
                    pen_type_str = 'condatTV', verbose_bl = verbose)
            else:
                # Use for 2D denoising
                denoise_cnf_obj = denoise.denoiseconf (log10_reg_lst = [-2, 2], nr_reg_int = 48, 
                    pen_type_str = 'TV', verbose_bl = verbose)
                adapt_adj = False
                try_denoise = False
                
            denoiser_obj = denoise.denoisepoisson (est_obj, denoise_cnf_obj)
            # Start the denoising
            denoiser_obj.denoise ()
            
            log_tune,valid_val = denoiser_obj.get_validation_loss()

            if accel and adapt_adj:
                index_min = np.argmin(valid_val)
                if verbose:
                    print('min index: %d'%index_min)
                # check if the optimiser chose end points
                # if so, adjust the tv values to test
                if accel_iter > max_iter:
                    # maximum number of iterations exceeded
                    try_denoise = False
                elif index_min == 0:
                    # reduce the range because the lowest tv penalty was chosen                
                    tv_reg[0] = tv_reg[0]*0.8
                    tv_reg[1] = tv_reg[0]*1.2
                    # check if the lowest range is below the tv lower limit
                    # if so, stop the loop
                    if tv_reg[0] < tv_lim[0]:
                        try_denoise=False
                elif index_min == len(log_tune)-1:
                    # increase the range because the greatest tv penalty was chosen                
                    tv_reg[0] = tv_reg[1]*0.8
                    tv_reg[1] = tv_reg[1]*1.2
                    # check if the range exceeds the predefined tv limits.
                    # if so, stop the loop here                
                    if tv_reg[1] > tv_lim[1]:
                        try_denoise=False
                
                else:
                    # a minimum was found that wasn't an end point.  Stop the loop
                    try_denoise = False
        
        ProfDenoise.profile[:,istart:iend] = denoiser_obj.getdenoised()/NumProf      
        
        
        tune_list.extend([[log_tune,valid_val]])  # store the results from the tuning parameters
        if plot_result:   
            plt.figure()
            plt.semilogy(ProfFit.flatten())
            plt.semilogy(denoiser_obj.getdenoised().flatten(),'--',linewidth=2)
            plt.title('Range Bins %d-%d'%(istart,iend))
            
            plt.figure()
            plt.plot(log_tune,valid_val)
            plt.title('Range Bins %d-%d'%(istart,iend))
            plt.show()
    
    
    if BG_Sub_Flag:
#        ProfDenoise.bg = ProfDenoise.bg
        ProfDenoise.profile = ProfDenoise.profile-ProfDenoise.bg[:,np.newaxis]
    
    ProfDenoise.ProcessingStatus.extend(['Applied time PTV denoising up to %.1f km'%(MaxAlt/1e3)])
    
    return ProfDenoise,tune_list        

def merge_hi_lo(hi_prof,lo_prof,lo_gain=0,plot_res=False):
    """
    merge_hi_lo(hi_prof,lo_prof,lo_gain=0,plot_res=False)    
    
    merge combined_hi and combined_lo gain profiles into a single profile.
    Nonlinear correction should be applied to the profiles prior to calling
    this function in order to capture nonlinear response as an uncertainty 
    in the hi gain channel.
    
    If no lo_gain adjustment is provided, the lo_gain is estimated

    hi_prof - high gain profile
    lo_prof - low gain profile    
    plot_res - if true, returns a scatter plot of the data before and after
        merging
    
        
    
    """
    
    if lo_gain == 0:
        # if lo_gain is not specified, estimate it using an optimizor
        errfun = lambda x: np.nansum((hi_prof.profile.flatten()-lo_prof.profile.flatten()*x)**2/(hi_prof.profile_variance.flatten()+lo_prof.profile_variance.flatten()*x**2))
        sol1D = scipy.optimize.minimize_scalar(errfun) #,disp=0,x0,maxfun=2000,eta=1e-5 , ,opt_iterations,opt_exit_mode
        lo_gain = sol1D['x']
        print('hi/lo combined gain estimate: %f'%lo_gain)
    lo_prof.gain_scale(lo_gain)
    
    
    
    combined = hi_prof.copy()
    # use low gain channel conditions
    combined_to_lo_prof = hi_prof.profile_variance > lo_prof.profile_variance  # variance of lo profile is lower
    combined_to_lo_prof = np.logical_or(combined_to_lo_prof,np.isnan(hi_prof.profile))  # OR combined_hi is nan
    combined_to_lo_prof = np.logical_and(combined_to_lo_prof,lo_prof.profile > 10)  # AND only use data point of combined_lo is significantly large
#    combined_to_lo = np.nonzero(np.logical_or( 
#        hi_prof.profile_variance > lo_prof.profile_variance,
#        np.isnan(hi_prof.profile)))
    combined_to_lo = np.nonzero(combined_to_lo_prof)
    combined.profile[combined_to_lo]= lo_prof.profile[combined_to_lo]
    combined.profile_variance[combined_to_lo]= lo_prof.profile_variance[combined_to_lo]
    combined.descript = 'Merged hi/lo gain combined channel'
    combined.label = 'Merged Combined Channel'
    
    if plot_res:
        plt.figure(); 
        plt.scatter(np.log10(lo_prof.profile.flatten()),np.log10(hi_prof.profile.flatten())); 
        plt.scatter(np.log10(lo_prof.profile.flatten()),np.log10(combined.profile.flatten())); 
        plt.plot([np.nanmin(np.log10(lo_prof.profile.flatten())),np.nanmax(np.log10(hi_prof.profile.flatten()))],[np.nanmin(np.log10(lo_prof.profile.flatten())),np.nanmax(np.log10(hi_prof.profile.flatten()))],'k--')
        plt.grid(b=True)
        plt.xlabel('Scaled Low Gain')
        plt.ylabel('High Gain/Merged')
        
    return combined,lo_gain
    
    
    
def AerosolBackscatter(MolProf,CombProf,CrossProf,Sonde,negfilter=True,eta_am=0.0,eta_ac=1.0,eta_mm=1.0,eta_mc=1.0,eta_x=0.0,gm=1.0):
    """
    Calculate the Aerosol Backscatter Coeffcient LidarProfiles: Molecular and Combined Channels
    Expects a 2d sonde profile that has the same dimensions as the Molecular and combined channels
    MolProf - molecular lidar profile (parallel polarized)
    CombProf - combined lidar profile (parallel polarized)
    CrossProf - cross polarized combined lidar profile
    Sonde - backscatter coefficient of molecular returns
    eta_am - transmission efficiency of aerosols into molecular channel (cross talk)
    eta_ac - transmission efficiency of aerosols into the combined channel
    eta_mm - transmission efficiency of molecules into the molecular channel (Rayleigh-Brillioun correction)
    eta_mc - transmission efficiency of molecules into the combined channel (Rayleigh-Brillioun correction)
    eta_x - parallel polarized cross talk into the cross polarized channel
    gm - molecular gain.  Set to 1.0 if adjusted prior to this function.
    Set negfilter = False to avoid filtering out negative values
    """
    
    ## Account for cross talk errors and Rayleigh-Brillioun efficiency effects
    # to obtain polarized molecular and aerosol backscatter signals
    dm = 0.00727345  # a priori known molecular depolarization
    mol = (MolProf*eta_ac-gm*CombProf*eta_am)/((dm-1)*gm*(eta_am*eta_mc-eta_ac*eta_mm))
    aer_par = (MolProf*eta_mc - gm*CombProf*eta_mm)/(gm*eta_am*eta_mc - gm*eta_ac*eta_mm)
    aer_per = (dm*MolProf*eta_ac*eta_mc+(CrossProf-CombProf*eta_x)*(gm*(eta_am*eta_mc-eta_ac*eta_mm))+dm*gm*(CrossProf*(-eta_am*eta_mc+eta_ac*eta_mm)+CombProf*eta_am*eta_mc*(-1+eta_x)-CombProf*eta_ac*eta_mm*eta_x))/((-1+dm)*gm*eta_ac*(-eta_am*eta_mc+eta_ac*eta_mm))
    
    comb_par = (MolProf*(eta_ac+(-1+dm)*eta_mc)-gm*CombProf*(eta_am+(-1+dm)*eta_mm))/((dm-1)*gm*(eta_am*eta_mc-eta_ac*eta_mm))
    
    mol.label = 'Molecular Signal'  
    mol.descript = 'Molecular signal estimated by inverting the cross coupling matrix'
    mol.profile_type = 'photon counts'
    
    aer_par.label = 'Parallel Aerosol Signal'  
    aer_par.descript = 'Parallel polarized aerosol signal estimated by inverting the cross coupling matrix'
    aer_par.profile_type = 'photon counts'
    
    aer_per.label = 'Perpendicular Aerosol Signal'  
    aer_per.descript = 'Perpendicular polarized aerosol signal estimated by inverting the cross coupling matrix'
    aer_per.profile_type = 'photon counts'
    
    comb_par.label = 'Combined Aerosol/Molecular Signal'  
    comb_par.descript = 'Parallel polarized combined aerosol and molecular signal estimated by inverting the cross coupling matrix'
    comb_par.profile_type = 'photon counts'
    
    param_profs= {'mol':mol,'aer_par':aer_par,'aer_per':aer_per,'comb_par':comb_par}
    
    Beta_AerBS = MolProf.copy()

    # calculate backscatter ratio
    BSR = (aer_par+aer_per+mol)/mol
    BSR.descript = 'Ratio of combined to molecular backscatter'
    BSR.label = 'Backscatter Ratio'
    BSR.profile_type = 'unitless'
    
#    Beta_AerBS.profile = (BSR-1)*beta_m_sonde[np.newaxis,:]    # only aerosol backscatter
    Beta_AerBS = (BSR-1)
#    Beta_AerBS.profile_variance = MolProf.profile_variance*(CombProf.profile)**2/(MolProf.profile)**4+CombProf.profile_variance*1/(MolProf.profile)**2  
    Beta_AerBS.multiply_prof(Sonde)  
    
    Beta_AerBS.descript = 'Calibrated Measurement of Aerosol Backscatter Coefficient in m^-1 sr^-1'
    Beta_AerBS.label = 'Aerosol Backscatter Coefficient'
    Beta_AerBS.profile_type = '$m^{-1}sr^{-1}$'
    
    dPart = aer_per/(aer_per+aer_par)
    dPart.descript = 'Propensity of Particles to depolarize (d).  This is not identical to the depolarization ratio.  See Gimmestad: 10.1364/AO.47.003795 or Hayman and Thayer: 10.1364/JOSAA.29.000400'
    dPart.label = 'Particle Depolarization'
    dPart.profile_type = 'unitless'
    
    
    if negfilter:
        Beta_AerBS.profile[np.nonzero(Beta_AerBS.profile <= 0)] = 1e-10;
    
    return Beta_AerBS,dPart,BSR,param_profs



