# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 08:50:06 2016

@author: mhayman
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import LidarProfileFunctions as lp
import scipy.interpolate

import glob


#filePath = '/scr/eldora1/MSU_h2o_data/2016/161209NF/18/'
#fileName1 = 'Online_Raw_Data.dat'  # molecular
#fileName2 = 'Offline_Raw_Data.dat'  # combined

use_mask = False
SNRmask = 3.0
countLim = 1.0



Day = 6
Month = 10
Year = 2017
HourLim = np.array([0,24])  # Limits on the processing time

tres = 1*60.0  # resolution in time points (0.5 sec)
zres = 10.0  # resolution in altitude points (7.5 m)

# index for where to treat the profile as background only
BGIndex = -100; # negative number provides an index from the end of the array
platform = 'ground' # 'ground' or 'airborne'.  If 'airborne' it needs an aircraft netcdf.
MaxAlt = 10e3


pol_xtalk = 0.015

kB = 1.3806504e-23;
c = 3e8

#basepath = '/scr/eldora1/HSRL_data/'  # old path - still works with link from HSRL_data to /hsrl/raw/
basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path

if Day < 10:
    DayStr = '0' + str(Day)
else:
    DayStr = str(Day)
    
if Month < 10:
    MonthStr = '0' + str(Month)
else:
    MonthStr = str(Month)

YearStr = str(Year)

    
FilePath0 = basepath + YearStr + '/' + MonthStr + '/' + DayStr + '/raw/'

SubFiles = glob.glob(FilePath0+'*.nc')
firstFile = True
for idir in range(len(SubFiles)):
    iTime = np.char.find(SubFiles[idir],'T')
    Hour = np.double(SubFiles[idir][iTime+1:iTime+3])
    if Hour >= np.floor(HourLim[0]) and Hour <= HourLim[1]:
        # Open the netcdf file and load the relevant variables
        f = netcdf.netcdf_file(SubFiles[idir], 'r')
        
        #TransEnergy = f.variables['total_energy'].data.copy();      # Transmitted pulse energy
        TransEnergy0 = lp.ncvar(f,'total_energy')
        
        # System for time array still needs work
        timeD = np.array(f.variables['DATA_time'].data);          # time array [year,month,day,hour,minute,second,msec,usec]
        timeDsec0 = timeD[:,3]*3600.0+timeD[:,4]*60.0+timeD[:,5]+timeD[:,6]*1e-3;    # convert time to seconds starting at the hour
        
        shot_count0 = f.variables['DATA_shot_count'].data.copy()     # number of laser shots in the time bin
        
        telescope_direction0 = f.variables['TelescopeDirection'].data.copy()  # indicates if telescope is pointed up or down (1 = up, -1? = down)
        telescope_locked0 = f.variables['TelescopeLocked'].data.copy()  # indicates if the telescope is locked - indicates possible tilted operation

        # Check if polarization (QWP Rotation) is in the netcdf.  Load if it is.
        if any('polarization' in s for s in f.variables):
            QWP = f.variables['polarization'].data.copy()               # QWP rotation angle
            # Calculate the QWP rotation rate to determine its status
            meanQWPrate = np.median(np.diff(QWP)/np.diff(timeDsec0)*180/np.pi)
            meanQWPvalue = np.median(QWP*180/np.pi)
            if meanQWPrate >= 10:
                QWP_Status = 'rotating'
            else:
                QWP_Status = 'fixed'
            print 'QWP Rotation Rate: %f deg/sec' %meanQWPrate
            print 'Average QWP Position: %f deg' %meanQWPvalue
        else:
            QWP_Status = 'fixed'
        
        print 'Processing %d UT' %Hour
        print 'Profile Integration Time: %f seconds' %np.median(np.diff(timeDsec0))
        print 'Processing QWP as %s' %QWP_Status
        
        cross_data = lp.profile_netcdf(f.variables['cross'])
        mol_data = lp.profile_netcdf(f.variables['molecular'])
        hi_data = lp.profile_netcdf(f.variables['combined_hi'])
        lo_data = lp.profile_netcdf(f.variables['combined_lo'])
#        print('%d,%d'(hi_data.data.shape[0],hi_data.data.shape[1]))
        f.close()
        
        # load profile data
        if firstFile:
            TransEnergy = TransEnergy0.copy()
            timeDsec = timeDsec0.copy()
            shot_count = shot_count0.copy()
            telescope_direction = telescope_direction0.copy()
            telescope_locked = telescope_locked0.copy()            
            
            CrossPol = lp.LidarProfile(cross_data,timeDsec0,label='Cross Polarization Channel',descript = 'Cross Polarization\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')
            RemCross = CrossPol.time_resample(delta_t=tres,update=True,remainder=True)
            
            Molecular = lp.LidarProfile(mol_data,timeDsec0,label='Molecular Backscatter Channel',descript = 'Parallel Polarization\nMolecular Backscatter Returns',bin0=47,lidar='GV-HSRL')
            RemMol = Molecular.time_resample(delta_t=tres,update=True,remainder=True)
            
            CombHi = lp.LidarProfile(hi_data,timeDsec0,label='High Gain Total Backscatter Channel',descript = 'Parallel Polarization\nHigh Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')
            RemCombHi = CombHi.time_resample(delta_t=tres,update=True,remainder=True)
            
            CombLo = lp.LidarProfile(lo_data,timeDsec0,label='Low Gain Total Backscatter Channel',descript = 'Parallel Polarization\nLow Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')            
            RemCombLo = CombLo.time_resample(delta_t=tres,update=True,remainder=True)
            
            
#            Molecular = lp.LidarProfile(mol_data.T,dt*np.arange(mol_data.shape[1])+Hour*3600,label='Molecular Backscatter Channel',descript = 'Unpolarization\nMolecular Backscatter Returns',bin0=0,lidar='DLB-HSRL')
#            RemMol = Molecular.time_resample(i=tres,update=True,remainder=True)
#            
#            CombHi = lp.LidarProfile(hi_data.T,dt*np.arange(hi_data.shape[1])+Hour*3600,label='Total Backscatter Channel',descript = 'Unpolarization\nHigh Gain\nCombined Aerosol and Molecular Returns',bin0=0,lidar='DLB-HSRL')
#            RemCom = CombHi.time_resample(i=tres,update=True,remainder=True)
            
            firstFile = False
            
            
        else:
            TransEnergy = np.concatenate((TransEnergy,TransEnergy0))
            timeDsec = np.concatenate((timeDsec,timeDsec0))
            shot_count = np.concatenate((shot_count,shot_count0))
            telescope_direction = np.concatenate((telescope_direction,telescope_direction0))
            telescope_locked = np.concatenate((telescope_locked,telescope_locked0))
            
            if np.size(RemMol.time) > 0:
                MolTmp = lp.LidarProfile(mol_data,timeDsec0,label='Molecular Backscatter Channel',descript = 'Parallel Polarization\nMolecular Backscatter Returns',bin0=47,lidar='GV-HSRL')
                MolTmp.cat_time(RemMol)
                RemMol = MolTmp.time_resample(delta_t=tres,update=True,remainder=True)
                Molecular.cat_time(MolTmp,front=False)
                
                CrossTmp = lp.LidarProfile(cross_data,timeDsec0,label='Cross Polarization Channel',descript = 'Cross Polarization\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')
                CrossTmp.cat_time(RemCross)
                RemCross = CrossTmp.time_resample(delta_t=tres,update=True,remainder=True)
                CrossPol.cat_time(CrossTmp,front=False)
                
                CombHiTmp = lp.LidarProfile(hi_data,timeDsec0,label='High Gain Total Backscatter Channel',descript = 'Parallel Polarization\nHigh Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')
                CombHiTmp.cat_time(RemCombHi)
                RemCombHi = CombHiTmp.time_resample(delta_t=tres,update=True,remainder=True)
                CombHi.cat_time(CombHiTmp,front=False)
                
                CombLoTmp = lp.LidarProfile(lo_data,timeDsec0,label='Low Gain Total Backscatter Channel',descript = 'Parallel Polarization\nLow Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')            
                CombLoTmp.cat_time(RemCombHi)
                RemCombLo = CombLoTmp.time_resample(delta_t=tres,update=True,remainder=True)
                CombLo.cat_time(CombLoTmp,front=False)
            else:
                MolTmp = lp.LidarProfile(mol_data,timeDsec0,label='Molecular Backscatter Channel',descript = 'Parallel Polarization\nMolecular Backscatter Returns',bin0=47,lidar='GV-HSRL')
                RemMol = MolTmp.time_resample(i=tres,update=True,remainder=True)
                Molecular.cat_time(MolTmp,front=False)
                
                CrossTmp = lp.LidarProfile(cross_data,timeDsec0,label='Cross Polarization Channel',descript = 'Cross Polarization\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')
                RemCross = CrossTmp.time_resample(delta_t=tres,update=True,remainder=True)
                CrossPol.cat_time(CrossTmp,front=False)
                
                CombHiTmp = lp.LidarProfile(hi_data,timeDsec0,label='High Gain Total Backscatter Channel',descript = 'Parallel Polarization\nHigh Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')
                RemCombHi = CombHiTmp.time_resample(delta_t=tres,update=True,remainder=True)
                CombHi.cat_time(CombHiTmp,front=False)
                
                CombLoTmp = lp.LidarProfile(lo_data,timeDsec0,label='Low Gain Total Backscatter Channel',descript = 'Parallel Polarization\nLow Gain\nCombined Aerosol and Molecular Returns',bin0=47,lidar='GV-HSRL')            
                RemCombLo = CombLoTmp.time_resample(delta_t=tres,update=True,remainder=True)
                CombLo.cat_time(CombLoTmp,front=False)
                



g1 = lp.load_geofile('/scr/eldora1/HSRL_data/2016/10/geofile_default_20161027T1610.geo')
#dg1 = lp.load_diff_geofile('/scr/eldora1/HSRL_data/2016/10/diff_default_geofile_20161027T1510.geo')
dg1 = lp.load_diff_geofile('/scr/eldora1/HSRL_data/2016/12/diff_geofile_20161219T1601.geo')
dgx = lp.load_diff_geofile('/scr/eldora1/HSRL_data/2016/10/cross_pol_diff_geofile_20161027T0000.geo',chan='cross')


#plt.figure(); 
#plt.plot(np.diff(QWP.data)/np.diff(timeDsec)*180/np.pi)


##### Pre-Processing of Profiles ####
EnergyNormFactor = 1.0/np.nanmean(TransEnergy.data)

CrossPol.slice_time(HourLim*3600)
CrossPol.nonlinear_correct(32e-9);
CrossPol.bg_subtract(BGIndex)
#CrossPol.energy_normalize(TransEnergy*EnergyNormFactor)
CrossPol.diff_geo_overlap_correct(dgx['cross']*dg1['hi'],geo_reference='mol')
CrossPol.geo_overlap_correct(1.0+g1)
CrossPol.range_correct()
CrossPol.slice_range(range_lim=[0,MaxAlt])
#CrossPol.conv(tres,zres)  # regrid by convolution
#CrossPol.time_resample(i=tres,update=True)  # regrid by binning
CrossPol.range_resample(delta_R=zres*7.5,update=True)

Molecular.slice_time(HourLim*3600)
Molecular.nonlinear_correct(28.4e-9);
Molecular.bg_subtract(BGIndex)
#Molecular.energy_normalize(TransEnergy*EnergyNormFactor)
Molecular.geo_overlap_correct(g1)
Molecular.range_correct();
Molecular.slice_range(range_lim=[0,MaxAlt])
#Molecular.conv(tres,zres)
#Molecular.time_resample(i=tres,update=True)  # regrid by binning
Molecular.range_resample(delta_R=zres*7.5,update=True)

CombHi.slice_time(HourLim*3600)
CombHi.nonlinear_correct(27.9e-9);
CombHi.bg_subtract(BGIndex)
#CombHi.energy_normalize(TransEnergy*EnergyNormFactor)
CombHi.diff_geo_overlap_correct(dg1['hi'],geo_reference='mol')
CombHi.geo_overlap_correct(g1)
CombHi.range_correct()
CombHi.slice_range(range_lim=[0,MaxAlt])
#CombHi.conv(tres,zres)
#CombHi.time_resample(i=tres,update=True)  # regrid by binning
CombHi.range_resample(delta_R=zres*7.5,update=True)

CombLo.slice_time(HourLim*3600)
CombLo.nonlinear_correct(21.3e-9);
CombLo.bg_subtract(BGIndex)
#CombLo.energy_normalize(TransEnergy*EnergyNormFactor)
CombLo.diff_geo_overlap_correct(dg1['lo'],geo_reference='mol')
CombLo.geo_overlap_correct(g1)
CombLo.range_correct()
CombLo.slice_range(range_lim=[0,MaxAlt])
#CombLo.conv(tres,zres)
#CombLo.time_resample(i=tres,update=True)  # regrid by binning
CombLo.range_resample(delta_R=zres*7.5,update=True)

#sigt = 3.0
#sigz = 1.0
#t = np.arange(-np.round(6*sigt),np.round(6*sigt))      
#z = np.arange(-np.round(6*sigz),np.round(6*sigz))  
#tt,zz = np.meshgrid(t,z)
#
#kconv = np.exp(-tt**2*1.0/(sigt**2)-zz**2*1.0/(sigz**2))
#kconv = kconv/(1.0*np.sum(kconv))
#Newprofile = scipy.signal.convolve2d(CombHi.profile,kconv,mode='same')
#Newprofile_variance = scipy.signal.convolve2d(CombHi.profile_variance,kconv**2,mode='same'


#plt.figure(); 
#plt.plot(np.sum(CombHi.profile,axis=0)/np.sum(CombLo.profile,axis=0),'.');
LoGain = np.mean((np.sum(CombHi.profile,axis=0)/np.sum(CombLo.profile,axis=0))[400:600])

#plt.figure(); 
#plt.plot(np.sum(CombHi.profile[3600:,:],axis=0)/np.sum(Molecular.profile[3600:,:],axis=0),'.');
#MolGain = 1.15*1.02;
#MolGain = 1.17  # after Dec 19, 2016
MolGain = 1.08 # Jan 9, 2017

CombLo.gain_scale(LoGain)
Molecular.gain_scale(MolGain)

lp.plotprofiles([CombHi,CombLo,Molecular,CrossPol])

### Grab Sonde Data
sondefilename = '/scr/eldora1/HSRL_data/2016/12/sondes.DNR.nc'
sonde_index = 2*Day
#(Man or SigT)
f = netcdf.netcdf_file(sondefilename, 'r')
TempDat = f.variables['tpSigT'].data.copy()  # Kelvin
PresDat = f.variables['prSigT'].data.copy()*100.0  # hPa - convert to Pa (or Man or SigT)
SondeTime = f.variables['relTime'].data.copy() # synoptic time: Seconds since (1970-1-1 00:00:0.0) 
SondeAlt = f.variables['htSigT'].data.copy()  # geopotential altitude in m
StatElev = f.variables['staElev'].data.copy()  # launch elevation in m
f.close()

TempDat[np.nonzero(np.logical_or(TempDat < 173.0, TempDat > 373.0))] = np.nan;
PresDat[np.nonzero(np.logical_or(PresDat < 1.0*100, PresDat > 1500.0*100))] = np.nan;

sonde_index = np.min([np.shape(SondeAlt)[0]-1,sonde_index])
# Obtain sonde data for backscatter coefficient estimation
Tsonde = np.interp(CombHi.range_array,SondeAlt[sonde_index,:]-StatElev[sonde_index],TempDat[sonde_index,:])
Psonde = np.interp(CombHi.range_array,SondeAlt[sonde_index,:]-StatElev[sonde_index],PresDat[sonde_index,:])

# note the operating wavelength of the lidar is 532 nm
beta_m_sonde = 5.45*(550.0/532.0)**4*1e-32*Psonde/(Tsonde*kB)

#plt.figure(); plt.semilogx(beta_m_sonde/np.nanmean(Molecular.profile,axis=0),Molecular.range_array)
Mol_Beta_Scale = 2.49e-11  # conversion from profile counts to backscatter cross section

BSR = (CombHi.profile+CrossPol.profile)/Molecular.profile
Depol = (CrossPol.profile-pol_xtalk*CombHi.profile)/(CombHi.profile+CrossPol.profile)

aer_depol = (-2*0.00365+BSR*Depol)/(BSR-1)
beta_bs = BSR*beta_m_sonde[np.newaxis,:]  # total backscatter including molecules
#aer_beta_bs = (BSR-1)*beta_m_sonde[np.newaxis,:]    # only aerosol backscatter
#aer_beta_bs[np.nonzero(aer_beta_bs <= 0)] = 1e-10;

aer_beta_gv = lp.Calc_AerosolBackscatter(Molecular,CombHi,Tsonde,Psonde)

extinction = -np.diff(np.log(Molecular.profile/beta_m_sonde[np.newaxis,:]),axis=1)*0.5

extinction2 = -np.diff(np.log(CombHi.profile/beta_bs),axis=1)*0.5

OptDepth = np.log(CombHi.profile/beta_bs)*0.5

#plt.figure(); 
#plt.imshow(np.log10(aer_beta_bs.T[-1:0:-1,:]));
##plt.imshow(np.log10(aer_beta_bs.T[-2000:0:-1,:]));
#plt.colorbar()
#plt.clim([-8,-3])

#mask_cmap = plt.cm.jet
#mask_cmap.set_bad('w',1.0)

if use_mask:
    CountMask = np.logical_or(Molecular.SNR() < countLim,CombHi.SNR() < countLim)
#    aer_beta_gv.profile = np.ma.array(aer_beta_gv.profile, mask=(np.logical_or(aer_beta_gv.SNR() < SNRmask,CountMask)))
    aer_beta_gv.profile = np.ma.array(aer_beta_gv.profile, mask=CountMask)
#    aer_beta_gv.profile[np.nonzero(aer_beta_gv.SNR() < SNRmask)] = np.nan

MskAer = np.ma.array(np.log10(aer_beta_gv.profile.T),mask=np.isnan(aer_beta_gv.profile.T))

plt.figure(figsize=(15,5)); 
#plt.imshow(np.log10(aer_beta_bs.T[-1:0:-1,:]));
plt.pcolor(aer_beta_gv.time/3600,aer_beta_gv.range_array*1e-3, np.log10(aer_beta_gv.profile.T));
plt.colorbar()
plt.clim([-9,-4])
plt.title(CombHi.lidar + ' Aerosol Backscatter Coefficient [$m^{-1}sr^{-1}$]')
plt.ylabel('Altitude [km]')
plt.xlabel('Time [UTC]')
plt.xlim(HourLim)

MolInt = Molecular.copy();
MolInt.time_integrate();
CombInt = CombHi.copy();
CombInt.time_integrate();
aer_beta_gv_int = lp.Calc_AerosolBackscatter(MolInt,CombInt,Tsonde,Psonde)
lp.plotprofiles([aer_beta_gv_int])

plt.figure(); 
plt.imshow(np.log10(aer_beta_gv.profile.T[-1:0:-1,:]));

#aer_beta_bs_gv = aer_beta_bs.copy()

plt.figure(); 
plt.imshow(aer_depol.T[-1:0:-1,:]);
#plt.imshow(aer_depol.T[-2000:0:-1,:]);
plt.colorbar()
plt.clim([0,0.7])