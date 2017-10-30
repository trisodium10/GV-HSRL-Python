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
    
    
def load_raw_data(start_time,stop_time,var_2d_list,var_1d_list,basepath = '/scr/eldora1/rsfdata/hsrl/raw/',verbose=True,as_prof=True):
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
    
    
    """
    
    var_1d_data = dict(zip(var_1d_list,[np.array([])]*len(var_1d_list)))
    var_2d_data = dict(zip(var_2d_list,[np.array([])]*len(var_2d_list)))
    timeD = np.array([])
    
    # build list of files falling on the requested dates
    day_iter =  start_time.date()
    SubFiles = []
    FileDate = []
    while day_iter <= stop_time.date():
        FilePath0 = basepath+day_iter.strftime('%Y/%m/%d/raw/')
        SubList0 = glob.glob(FilePath0+'*data*.nc')
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
        if File_Time >= start_time and File_Time <= stop_time:
            
            f = nc4.Dataset(SubFiles[idir],'r')
            
            for var in var_1d_data.keys():
                if any(var in s for s in f.variables):
                    var_data = lp.ncvar(f,var)
                    var_1d_data[var] = np.concatenate((var_1d_data[var],var_data)) 
                
            # System for time array still needs work
            timeD0 = np.array(f.variables['DATA_time'][:]).astype(np.int)
            timeD0[:,6] = timeD0[:,6]*1e3+timeD0[:,7]
            time_dt0 = [datetime.datetime(*x) for x in timeD0[:,:7]]
            time_sec0 = np.array([(x-time_dt0[0]).total_seconds() for x in time_dt0])
            
            if len(timeD)>0:
                timeD = np.vstack((timeD,np.array(f.variables['DATA_time'][:].astype(np.int))));          # time array [year,month,day,hour,minute,second,msec,usec]
            else:
                timeD = np.array(f.variables['DATA_time'][:]).astype(np.int)
            
            if verbose:
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
                bin0=47,lidar='GV-HSRL',StartDate=start_time.date())]
            prof_names = prof_names + [var]
            
        profs_2d = dict(zip(prof_names,prof_list))
    else:
        profs_2d = var_2d_data
    
    return time_list, var_1d_data, profs_2d