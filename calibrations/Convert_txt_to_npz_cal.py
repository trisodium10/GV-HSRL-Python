# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 19:29:09 2017

@author: mhayman

Converts text based cal file into an npz for running the NCAR
GV-HSRL processing scripts

"""

import numpy as np
import datetime

#fname = 'diff_default_geofile_20171025T1805.geo'
fname  = 'diff_default_geofile_20120201T0000.geo'

#read_path = '/Users/mhayman/Documents/Python/Lidar/GV-HSRL-Python/calibrations/cal_files/'
read_path = '/scr/eldora1/HSRL_data/2015/07/'

write_path = '/Users/mhayman/Documents/Python/Lidar/' #GV-HSRL-Python/calibrations/cal_files/'



cal = np.loadtxt(read_path + fname)

if '.geo' in fname:
    f= open(read_path+fname, "r")
    file_list = f.readlines()[:2]
    f.close()
    date_str = file_list[0].split('profile')[1]
    effective_date = datetime.datetime.strptime(fname[-17:-4],'%Y%m%dT%H%M')
    cal_date = datetime.datetime.strptime(date_str.split(' -->')[0],' data from %d-%b-%y %H:%M')
    stop_time = datetime.datetime.strptime(cal_date.strftime('%d-%b-%y ')+date_str.split(' -->')[1],'%d-%b-%y %H:%M UTC\n')
    save_file_name = write_path+'diff_geo_GVHSRL'+ effective_date.strftime('%Y%m%d')
    np.savez(save_file_name,hi_diff_geo=cal[:,1],lo_diff_geo=cal[:,2],start_time = cal_date,stop_time=stop_time)
    