# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 09:09:53 2018

@author: mhayman
"""

import os

#/scr/rain1/rsfdata/projects/socrates/hsrl/raw/2018/01/02/raw/gvhsrl_201801

#basepath = '/scr/eldora1/HSRL_data/'  # old path - still works with link from HSRL_data to /hsrl/raw/
#basepath = '/scr/eldora1/rsfdata/hsrl/raw/'  # new absolute path

basepath = '/Users/mhayman/Documents/HSRL/GVHSRL_data/'  # local computer data path
#basepath = '/scr/rain1/rsfdata/projects/socrates/hsrl/raw/' #SOCRATES basepath

# Server path (FL)
#aircraft_basepath = {
#    'CSET':'/scr/raf_data/CSET/24_Mar_17_BU/',
#    'SOCRATES':'/scr/raf_data/SOCRATES/' #SOCRATEStf01.nc       
#    } 
    
# Local Paths
aircraft_basepath = {
    'CSET':'/Users/mhayman/Documents/HSRL/aircraft_data/',
    'SOCRATES':'/Users/mhayman/Documents/HSRL/aircraft_data/' #SOCRATEStf01.nc       
    } 

save_data_path = os.environ['HOME'] + '/Documents/HSRL/GVHSRL/data/'
save_plots_path = os.environ['HOME'] + '/Documents/HSRL/GVHSRL/plots/'

software_path = os.environ['HOME']+'/Documents/Python/Lidar/GV-HSRL-Python/'

#Processor = os.environ['HOME']+'/Documents/Python/Lidar/GV-HSRL-Python/processors/Airborne_Processor_GVHSRL.py'
