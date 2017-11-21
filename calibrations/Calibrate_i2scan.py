# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 09:41:02 2017
@author: mhayman
"""





##read theoretical i2 transmission file
#    try:
#         [header,i2t]=ru.readascii(locate_file('dfb_i2_spectra.txt'))
#         dfb_to_i2_scale = rs_constants['i2_absorption_scale_factor']
#         i2t_trans=10**(dfb_to_i2_scale*i2t[:,1])
#         i2t_freq = i2t[:,0]*1e9         
#    except IOError: #file either not found by locate, or error occurred 
#         [header,i2t]=ru.readascii(locate_file('I2cell_272_31_extended_freq.txt'))
#         i2t_freq  = i2t[:,0]*1e9
#         i2t_trans = i2t[:,2]
#
#
#    #normalize theoretical spectrum to 1--ie eliminate continuum absorption
#    lo_freq_lmt = 2e9
#    hi_freq_lmt = 4e9
#    mask = np.array([ ((f >= lo_freq_lmt) and (f < hi_freq_lmt)) for f in i2t_freq])
#    print 'normalizing i2 transmission between ', lo_freq_lmt/1e9, ' and ',hi_freq_lmt/1e9, ' GHz'
#    norm = nanmean(i2t_trans[mask])
#    i2t_trans = i2t_trans/norm

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
#from scipy import optimize
import netCDF4 as nc4

import glob

import os

import datetime

def Tetalon(x,nu):
    T = (1-x[0])**2/((1-x[0])**2+4*x[0]*np.sin(x[2]*freq+x[1])**2)
    return T

def Fit_Error(x,comb,mol,freq,i2,i2f,plotres=False):
    """
    x[0] frequency multiplier
    x[1] frequency offset
    x[2] etalon reflection ceoff
    x[3] etalon frequency offset
    x[4] etalon phase multiplier
    x[5] molecular multiplier
    x[6] combined multiplier
    """
    freq2 = x[0]*freq+x[1]
#    Cfft = np.fft.fft(np.log(comb))
#    Mfft = np.fft.fft(np.log(mol))
    Te = Tetalon(x[2:5],freq2)
    i2spec = np.interp(freq2,i2f,i2)
    Mideal = x[5]*i2spec*Te
    Cideal = x[6]*Te
#    Mideal_fft = np.fft.fft(np.log(Mideal))
#    Cideal_fft = np.fft.fft(np.log(Cideal))
    
    if plotres:
#        plt.figure()
#        plt.plot(np.real(Mfft))
#        plt.plot(np.real(Mideal_fft))
#        
#        plt.figure()
#        plt.plot(np.imag(Mfft))
#        plt.plot(np.imag(Mideal_fft))
#        
#        plt.figure()
#        plt.plot(np.real(Cfft))
#        plt.plot(np.real(Cideal_fft))
#        
#        plt.figure()
#        plt.plot(np.imag(Cfft))
#        plt.plot(np.imag(Cideal_fft))
        
        plt.figure()
        plt.semilogy(freq2,mol)
        plt.semilogy(freq2,Mideal)
        
        plt.figure()
        plt.semilogy(freq2,comb)
        plt.semilogy(freq2,Cideal)
    
#    Error = np.nansum((np.log(comb)-np.log(Cideal))**2)
#    Error = np.nansum((Cfft-Cideal_fft)*((Cfft-Cideal_fft).conj()))
#    Error = np.nansum((np.real(Cfft)-np.real(Cideal_fft))**2)+np.nansum((np.imag(Cfft)-np.imag(Cideal_fft))**2)
#    Error = -1*np.nansum(Mfft*Mideal_fft.conj())-np.nansum(Cfft*Cideal_fft.conj())
#    Error = np.real(np.nansum((Mfft-Mideal_fft)*(Mfft-Mideal_fft).conj()+(Cfft-Cideal_fft)*(Cfft-Cideal_fft).conj()))
    Error = np.nansum((np.log(mol)-np.log(Mideal))**2)+np.nansum((np.log(comb)-np.log(Cideal))**2)
    return Error
    
def get_i2_write_values(x,comb,mol,freq,i2,i2f,plotres=False):
    """
    x[0] frequency multiplier
    x[1] frequency offset
    x[2] etalon reflection ceoff
    x[3] etalon frequency offset
    x[4] etalon phase multiplier
    x[5] molecular multiplier
    x[6] combined multiplier
    """
    freq2 = x[0]*freq+x[1]
    Te = Tetalon(x[2:5],freq2)
    molecular = x[5]/x[6]*np.interp(i2f,freq2,mol)
    combined = np.interp(i2f,freq2,comb)
    i2_measured = np.interp(i2f,freq2,mol/Te)
    i2_theory = i2
    Te2 = Tetalon(x[2:5],i2f)
    if plotres:
        plt.figure()
        plt.plot(i2f,combined)
        plt.plot(i2f,molecular)
        plt.plot(i2f,i2_measured)
        plt.plot(i2f,i2_theory)
    return np.vstack((i2f,combined,molecular,i2_measured,i2_theory)).T,Te2

# input and raw_input are the same for python 3 and 2 respectively
# this makes it so input always accepts a string
try:
    input=raw_input
except NameError:
    pass


print("Enter I2 Calibration:")
year_in = np.int(input("Year: "))
month_in = np.int(input("Month (#): "))
day_in = np.int(input("Day: "))

cal_date = datetime.datetime(year_in,month_in,day_in)

basepath = '/scr/eldora1/HSRL_data/'

FilePath0 = basepath+cal_date.strftime('%Y/%m/%d/raw/')

file_list = glob.glob(FilePath0+'*calibration*.nc')
print('Searching: \n'+FilePath0+'\n')
print('Found Calibration Files:')
for ai in range(len(file_list)):
    print('%d.)  '%ai+file_list[ai].split(FilePath0,1)[1]+'\n')

file_select = np.int(input('Select File Number (-1==quit): '))

if file_select != -1 and file_select < len(file_list):
    


    """
    int LIDARCTL_mode(time) ;
    		LIDARCTL_mode:long_name = "Lidar Mode Bits" ;
    		LIDARCTL_mode:bit_9 = "Iodine Seek" ;
    		LIDARCTL_mode:bit_10 = "Energy Monitor ND Filter In" ;
    		LIDARCTL_mode:bit_11 = "Cross channel blocked" ;
    		LIDARCTL_mode:bit_0 = "Iodine Scan" ;
    		LIDARCTL_mode:bit_1 = "Etalon Scan" ;
    		LIDARCTL_mode:bit_2 = "Iodine Lock" ;
    		LIDARCTL_mode:bit_3 = "Etalon Lock" ;
    		LIDARCTL_mode:bit_4 = "part of a Calibration Scan" ;
    		LIDARCTL_mode:bit_5 = "Transmitter or Receiver Attenuated for Calibration" ;
    		LIDARCTL_mode:bit_5_notes = "AHSRL transmission attenuator is low or GVHSRL (and others) receiver attenuating filter is in.  Can be verified which by checking the transmitted energy monitor" ;
    		LIDARCTL_mode:bit_6 = "Attenuator is moving" ;
    		LIDARCTL_mode:bit_7 = "Neutral Density Filter In for wide scan" ;
    		LIDARCTL_mode:bit_8 = "Part of an Etalon-tuning Iodine Scan" ;
    		LIDARCTL_mode:canonical_device_name = "lidarmode" ;
    
    """
    
    
    
        
    
    
    phase_to_GHz = -1/5.0
    
    tstart = 1.5e-6
    tend = 2.1e-6
    
    i_cal_0 = np.round(tstart/50e-9).astype(np.int)
    i_cal_1 = np.round(tend/50e-9).astype(np.int)
    
    print('Calibrating using '+file_list[file_select])
    f = nc4.Dataset(file_list[file_select],'r')
    #f = netcdf.netcdf_file('/scr/eldora1/HSRL_data/2017/10/18/raw/gvhsrl_20171018T190138_calibration_fl1.nc', 'r')
#    f = nc4.Dataset('/scr/eldora1/HSRL_data/2017/10/18/raw/gvhsrl_20171018T190138_calibration_fl1.nc','r')
    
    comb = f.variables['combined_hi'][:].copy()
    mol = f.variables['molecular'][:].copy()
    
    lidarctl = f.variables['LIDARCTL_mode'][:].copy()
    
    scan_start = f.variables['DATA_first_time'][0].copy()
    scan_end = f.variables['DATA_last_time'][-1].copy()
    
    interferometer = f.variables['interferometer'][:].copy()
    
    scan_start_dt = datetime.datetime(scan_start[0],scan_start[1],scan_start[2],scan_start[3],scan_start[4])
    scan_stop_dt = datetime.datetime(scan_end[0],scan_end[1],scan_end[2],scan_end[3],scan_end[4])
    
    #l3volt = f.variables['l3cavityvoltage'][:].copy()
    #seed_temp = f.variables['DATA_seed_laser_temp'][:].copy()
    #seed_hist = f.variables['seedcontroller_histogram'][:].copy()
    f.close()
    
    molcal = np.sum(mol[:,i_cal_0:i_cal_1],axis=1)  # cal-pulse in the molecular channel
    combcal = np.sum(comb[:,i_cal_0:i_cal_1],axis=1)
    
    i2scan = lidarctl&1 == 1            # is the system doing a scan?
    i2wide = lidarctl&(1|128) == (1|128)  # is the ND filter in for the scan?
    i2narrow = lidarctl&(1|128) == 1  # is the ND filter in for the scan?
    
    i_wide = np.nonzero(i2wide)[0]
    i_narrow = np.nonzero(i2narrow)[0]
    
    #plt.figure();
    #plt.semilogy(molcal[i_wide])
    #plt.semilogy(molcal[i_narrow]*1e-3)
    
    #  Load I2 spectrum
    i2file = os.path.abspath(__file__+'/../../reference_files/')+'/I2cell_272_31_extended_freq.txt'
#    i2file = '/h/eol/mhayman/PythonScripts/HSRL_Processing/NewHSRLPython/I2cell_272_31_extended_freq.txt'
    i2spec = np.loadtxt(i2file)  #,skiprows=4
    
    save_ez_file_path = os.path.abspath(__file__+'/../cal_files/')+'/'
#    save_file_path = '/h/eol/mhayman/HSRL/hsrl_processing/hsrl_configuration/projDir/calfiles/'
    
    MZfft = np.fft.fft(interferometer,axis=1)
    imin = 15
    imax = 35
    ipk = np.argmax(np.abs(MZfft[:,imin:imax]),axis=1)+imin
    MZphase = np.angle(MZfft)[np.arange(MZfft.shape[0]),ipk]
    
    MZphase = np.unwrap(MZphase)
    phase_wide = MZphase[i_wide]
    phase_narrow = MZphase[i_narrow]
    
    #phase_wide = np.unwrap(MZphase[i_wide])
    #phase_wide = phase_wide-phase_wide[np.int(phase_wide.size/2.0)]
    #phase_narrow = np.unwrap(MZphase[i_narrow])
    #phase_narrow = phase_narrow - phase_narrow[np.int(phase_narrow.size/2.0)]
    
    #plt.figure();
    #plt.semilogy(phase_wide,molcal[i_wide])
    #plt.semilogy(phase_narrow,molcal[i_narrow]*1e-3)
    
    dphase = np.mean(np.diff(phase_narrow))
    phase_wide_interp = np.arange(phase_wide[0],phase_wide[-1],dphase)
    molwide = np.interp(phase_wide_interp,phase_wide,molcal[i_wide])
    norm_factor = np.max(molwide)
    
    ip0 = np.argmin(np.abs(phase_wide_interp-phase_narrow[0]))
    if phase_wide_interp[ip0] < phase_narrow[0]:
        ip0 = ip0+1
    ip1 = np.argmin(np.abs(phase_wide_interp-phase_narrow[-1]))
    if phase_wide_interp[ip1] > phase_narrow[-1]:
        ip1 = ip1-1
    phase_nar_interp = np.arange(phase_wide_interp[ip0],phase_wide_interp[ip1],dphase)
    #phase_nar_interp = np.arange(phase_narrow[0],phase_narrow[-1],dphase)
    molnar = np.interp(phase_nar_interp,phase_narrow,molcal[i_narrow])
    
    corMolNar = np.log(molnar)
    corMolNar[np.isinf(corMolNar)] = -5.0
    corMolWide = np.log(molwide)
    corMolWide[np.isinf(corMolWide)] = -5.0
    
    cormol = np.correlate(corMolNar,corMolWide)
    imid = np.int(cormol.size/2)
    ipkcor = -30+np.argmax(cormol[imid-30:imid+30])
    
    #plt.figure(); 
    #plt.semilogy(phase_nar_interp+ipkcor*dphase,molnar*1e-3/norm_factor)
    #plt.semilogy(phase_wide_interp,molwide/norm_factor)
    
    iswap = np.nonzero(molnar*1e-3>molwide[ip0:ip1+1])[0]
    
    mol_merge = molwide
    mol_merge[iswap[0]+ip0:iswap[-1]+ip0] = 1e-3*molnar[iswap[0]:iswap[-1]]
    mol_merge = mol_merge/norm_factor
    
    #plt.figure()
    #plt.semilogy(phase_wide_interp,mol_merge)
    
    
    
    freq = phase_wide_interp[::-1]*phase_to_GHz
    mol_merge = mol_merge[::-1]
    combwide = np.interp(phase_wide_interp,phase_wide,combcal[i_wide])
    comb = combwide[::-1]
    comb = comb/np.max(comb)
    
    #plt.figure()
    #plt.semilogy(i2spec[:,0],i2spec[:,2])
    #plt.semilogy(freq,mol_merge)
    
    
    """
        x[0] frequency multiplier
        x[1] frequency offset
        x[2] etalon reflection ceoff
        x[3] etalon frequency offset
        x[4] etalon phase multiplier
        x[5] molecular multiplier
        x[6] combined multiplier
    """
    
    
    
    #Tetalon = (1-R)**2/((1-R)**2+4*R*np.sin(delta))
    #Fit_Error(x,comb,mol,freq,i2,i2f,plotres=False):
    fit_fun = lambda x: Fit_Error(x,comb,mol_merge,freq,i2spec[:,2],i2spec[:,0])
    x0 = np.array([0.85,-0.8,0.90,0.0,0.01,1.0,1.0])
    bnds = np.zeros((x0.size,2))
    bnds[0,1] = 10
    bnds[1,0] = -10
    bnds[1,1] = 10
    bnds[2,1] = 1
    bnds[3,0] = -10
    bnds[3,1] = 10
    bnds[4,1] = 10
    bnds[5,1] = 100
    bnds[6,1] = 100
    
    sol = scipy.optimize.fmin_slsqp(fit_fun,x0,bounds=bnds)
    Fit_Error(sol,comb,mol_merge,freq,i2spec[:,2],i2spec[:,0],plotres=True)
    
    
    print('Phase to frequency conversion: %f'%(phase_to_GHz*sol[0]))
    print('Phase to frequnecy offset: %f'%sol[1])
    
    write_data,Te = get_i2_write_values(sol,comb,mol_merge,freq,i2spec[:,2],i2spec[:,0],plotres=True)
    
    plt.show(block=False)
    
    save_i2_str = input("Save i2 cal? [Y/N]  ")
    if save_i2_str[0] == 'Y' or save_i2_str[0] == 'y':  
        save_ez_filename =  'i2-default-scan-'+scan_start_dt.strftime('%Y%m%dT%H%M')
        np.savez(save_ez_file_path+save_ez_filename,freq_mult = sol[0], \
            freq_offset = sol[1], etalon_refl = sol[2], etalon_freq_offset = sol[3], \
            etalon_phase_mult = sol[4], molecular_mult = sol[5], combined_mult = sol[6], \
            freq=write_data[:,0],combined_scan=write_data[:,1],mol_scan=write_data[:,2],\
            i2_theory=write_data[:,4],Tetalon = Te,created_str = datetime.datetime.today().strftime('%d-%b-%y at %H:%m UT'))
        
        
        header_str = ''
        header_str = header_str+'Calibration scan as function of frequency offset from I2 line\n'
        header_str = header_str+'calibration scan data aquired on ' + scan_start_dt.strftime('%d-%b-%y at %H:%M UT') + '\n'
        header_str = header_str+'file created on '  + datetime.datetime.today().strftime('%d-%b-%y at %H:%m UT') + '\n'
        header_str = header_str+'t_begin_cal_pulse= %1.2e ;  start time of cal pulse (sec).\n'%tstart
        header_str = header_str+'t_end_cal_pulse=   %1.2e ;  end time of cal pulse (sec).\n'%tend
        header_str = header_str+'pulse_durration=   5.00e-08 ;  laser pulse durration (sec).\n'
        header_str = header_str+'ratio of mol to combined channel gain =   %.1f \n' %(sol[5]/sol[6])
        header_str = header_str+'Cam = %1.3e \n'%(np.min(mol_merge)*sol[6]/sol[5])
        header_str = header_str+'Min iodine transmission =    %1.1e,  1/(min_trans) =  %1.1e\n\n'%(np.min(mol_merge),1/np.min(mol_merge))
        header_str = header_str+'freq(GHz)  combined  molecular i2_measured i2_theory'
        
#        # may need 'i2-default-scan' in name
#        save_cal_file = 'i2-default-scan-'+scan_start_dt.strftime('%Y%m%dT%H%M')+'.cal'
##        save_file_path = '/h/eol/mhayman/HSRL/hsrl_processing/'
#        save_file_path = '/h/eol/mhayman/HSRL/hsrl_processing/hsrl_configuration/projDir/calfiles/'
#        np.savetxt(save_file_path+save_cal_file,write_data,fmt='%.5f    ',header=header_str)

else:
    print('Unrecognized file index')
    print('Terminating Program')

#i2calfile = '/h/eol/mhayman/HSRL/hsrl_processing/hsrl_configuration/projDir/calfiles/i2-default-scan-20150601T0000.cal'
#i2cal = np.loadtxt(i2calfile)

#plt.figure()
#plt.plot(i2cal[:,0],i2cal[:,1])
#plt.plot(i2cal[:,0],i2cal[:,2])
#plt.plot(i2cal[:,0],i2cal[:,3])
#plt.plot(i2cal[:,0],i2cal[:,4])