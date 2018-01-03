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
    
    
def load_raw_data(start_time,stop_time,var_2d_list,var_1d_list,basepath = '/scr/eldora1/rsfdata/hsrl/raw/',verbose=True,as_prof=True,loadQWP='fixed'):
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
    
    
    """
    
    var_1d_data = dict(zip(var_1d_list,[np.array([])]*len(var_1d_list)))
    var_2d_data = dict(zip(var_2d_list,[np.array([])]*len(var_2d_list)))
    timeD = np.array([])
    
    start_time_hr = datetime.datetime(year=start_time.year,month=start_time.month,day=start_time.day,hour=start_time.hour)    
    
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
        iTime = np.char.find(SubFiles[idir],'T')
        Hour = np.double(SubFiles[idir][iTime+1:iTime+3])
        File_Time = datetime.datetime(*FileDate[idir].timetuple()[:4])+datetime.timedelta(hours=Hour)
        # check if this file is in the search time period
        if File_Time >= start_time_hr and File_Time <= stop_time:
            
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
                    print( 'Processing %d UT' %Hour)
                    print( '['+SubFiles[idir]+']')
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
    #timeD = timeD[:,:8]
    timeD[:,6] = timeD[:,6]*1e3+timeD[:,7]
    time_dt = [datetime.datetime(*x) for x in timeD[:,:7]]
    time_dt0 = datetime.datetime(time_dt[0].year,time_dt[0].month,time_dt[0].day)
    time_sec = np.array([(x-time_dt0).total_seconds() for x in time_dt])
    
    time_dt = np.array(time_dt)
    
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
    
    return time_list, var_1d_data, profs_2d

def load_aircraft_data(filename,var_list):
    """
    Loads data from aircraft netcdf
    filename - string with path and filename to load
    var_list - strings of the netcdf variables to load
    """    

    var_data = dict(zip(var_list,[np.array([])]*len(var_list)))
    f = nc4.Dataset(filename,'r')
            
    for var in var_data.keys():
        if any(var in s for s in f.variables):
            data = lp.ncvar(f,var)
            if len(var_data[var]) > 0:
                var_data[var] = np.concatenate((var_data[var],data)) 
            else:
                var_data[var] = data.copy()
                
    print('Aircraft Time Data: %s' %f.variables['Time'].units)
    f.close()
    
    return var_data

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
        iremain = np.int(np.max(itime))  # the remainder starts at the maximum bin where data exists
#                if not remainder and iremain < tedges.size:
#                    iremain = iremain+1  # make sure to include the last bin of data if remainder isn't being used.
        iremainList = np.nonzero(var_time > tedges[iremain-1])[0]
        iprofstart = np.int(np.max(np.array([1,np.min(itime)])))
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
            
def get_TP_from_aircraft(air_data,profile):
    """
    Returns a temperature and pressure profile corresponding to
    dimensions of profile
    The profile should already be converted to altitude coordinates from range
    
    The aircraft temperature and pressure measurements are contained
    in air_data
    """
    
    # make sure the time axes of the aircraft data aline with the 
    # lidar profile data
    air_data_int = interp_aircraft_data(profile.time,air_data)
    aircraft_temp = air_data_int['ATX']+273.15  # convert aircraft temperature to K 
    b_T = (aircraft_temp + air_data_int['GGALT']*0.0065)[:,np.newaxis]
    TempAir = b_T-0.0065*profile.range_array[np.newaxis,:]
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