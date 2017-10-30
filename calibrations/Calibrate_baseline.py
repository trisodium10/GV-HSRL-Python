# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 12:47:00 2017

@author: mhayman
"""

import sys
import os

# add the path to GVHSRLlib manually
library_path = os.path.abspath(__file__+'/../../libraries/')
print(library_path)
if library_path not in sys.path:
    sys.path.append(library_path)

import numpy as np
import datetime as datetime

import GVHSRLlib as gv

"""
Write a baseline file with all zeros
"""



save_path_ez = os.path.abspath(__file__+'/../cal_files/')+'/'

#save_file_path = '/Users/mhayman/Documents/Python/Lidar/'
#save_file_path = '/h/eol/mhayman/HSRL/hsrl_processing/hsrl_configuration/projDir/calfiles/'
save_file_path = save_path_ez


sg_win = 11
sg_order = 5
        
year_in = 2017
month_in = 10
day_in = 19
start_hr = 20
stop_hr = 2

print('This program only writes empty (zeros) baseline files')

print('Default Test Date:')
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


var_1d_list = ['total_energy','RemoveLongI2Cell'\
    ,'TelescopeDirection','TelescopeLocked','polarization']  # 'DATA_shot_count'


var_2d_list = ['cross','molecular','combined_hi','combined_lo']



save_file_name = 'baseline_correction_'+cal_start.strftime('%Y%m%dT%H%M')+'.blc'


write_data = np.zeros((4000,7))
write_data[:,0] = np.arange(write_data.shape[0])

# write out zeros for baseline correction
#file = open(save_file_path+save_file_name,'w') 
#file.write('#profile data from ' + cal_start.strftime('%d-%b-%y %H:%M')+' -->'+cal_stop.strftime('%H:%M UTC'))
#file.write('#ave energy per shot= 0.06000    mJ')
#header_str = 'bin_num\tcombined_hi\tcombined_lo\tmolecular\tcrosspol\tmol_I2A\tcomb_1064'
#file.close()

header_str = 'profile data from ' + cal_start.strftime('%d-%b-%y %H:%M')+' -->'+cal_stop.strftime('%H:%M UTC') + '\n'
header_str = header_str + 'ave energy per shot= 0.06000    mJ\n'
header_str = header_str + 'bin_num\tcombined_hi\tcombined_lo\tmolecular\tcrosspol\tmol_I2A\t\tcomb_1064'

np.savetxt(save_file_path+save_file_name,write_data,fmt='\t%d\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e',header=header_str,comments='#')
