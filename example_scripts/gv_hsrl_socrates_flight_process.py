# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 07:59:40 2018

@author: mhayman

Runs two python scripts.
The first processes flight data and saves it as a netcdf.
The second plots the netcdf data and saves the plots.

These scripts are configured for SOCRATES processing.  The user will
be prompted to select which SOCRATES flight to process.


"""
import os

ProcessFile = os.path.abspath(__file__+'/../')+'/gv_hsrl_socrates_process_hour.py'
PlotFile = os.path.abspath(__file__+'/../')+'/gv_hsrl_socrates_plot_ncdata.py'
PlotFile2 = os.path.abspath(__file__+'/../')+'/gv_hsrl_socrates_plot_ncdata_fiveminute.py'

# Process the data (and save it)
exec(open(ProcessFile).read())

# Plot the data (and save the plots) in 1 hour rolling profiles
exec(open(PlotFile).read())

# Plot the data (and save the plots) in five minute profiles with kft altitude
exec(open(PlotFile2).read())

