#!/usr/bin/env python

import numpy as np
import fileinput
from decimal import *
import datetime
from math import ceil,floor,cos,sin,acos,tan,exp,fabs,sqrt
import shlex
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
import os
import shutil
import sys
from scipy.interpolate import interp1d
import time
from matplotlib.dates import MonthLocator, YearLocator, DateFormatter, drange
from operator import add, sub

from netCDF4 import Dataset, date2num

import calendar
import datetime as dt
from dateutil.relativedelta import relativedelta

##########################################################################################################
# add sites to this array for model runs
#sites=["BAN", "FNS", "K34", "K67", "K77", "K83", "PDG", "RJA", "CAX04", "CAX06", "KEN01", "KEN02", "TAM05", "TAM06", "NVX"]
sites=["BAN","FNS","CAX","K34","K67","K77","K83","PDG","RJA"] # EDIT THIS, add site names to this array

# Can be lba or gem
obs_data="lba"

# site: [lat, lon]
sitelocs = {
"CAX04": [-1.716, -51.457],
"CAX06": [-1.737, -51.462],
"KEN01": [-16.016, -62.730],
"KEN02": [-16.016, -62.730],
"TAM05": [-12.831, -69.271],
"TAM06": [-12.839, -69.296],
"NVX": [-14.70, -52.35],
"BAN": [-9.82, -50.15],
"FNS": [-10.76, -62.36],
"K34": [-2.61, -60.21],
"K67": [-2.85, -54.97],
"K77": [-3.01, -54.54],
"K83": [-3.05, -54.93],
"PDG": [-21.62, -47.63],
"RJA": [-10.08, -62.36],
"CAX": [-1.748, -51.454]}

# Note on CAX - lat lon taken from the site info from the LCB dataset:
# file:///exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/FLUX_LBA_ECO/CD32_BRAZIL_FLUX_NETWORK_1174/guide/CD32_Brazil_Flux_Network.html

# start date and number of daily time steps
site_dates = {
"BAN": [dt.date(2004, 1, 2), 1033],
"FNS": [dt.date(1999, 1, 2), 1094],
"K34": [dt.date(2003, 1, 2), 1017],
"K67": [dt.date(2002, 1, 2), 685],
"K77": [dt.date(2001, 1, 2), 1824],
"K83": [dt.date(2001, 1, 2), 952],
"PDG": [dt.date(2002, 1, 2), 728],
"RJA": [dt.date(2000, 2, 3), 953],
"CAX": [dt.date(1999, 1 , 1), 731],
"CAX04": [dt.date(2005, 1, 1), 4383],
"CAX06": [dt.date(2005, 1, 1), 4383],
"KEN01": [dt.date(2005, 1, 1), 4383],
"KEN02": [dt.date(2005, 1, 1), 4383],
"TAM05": [dt.date(2005, 1, 1), 4383],
"TAM06": [dt.date(2005, 1, 1), 4383],
"NVX": [dt.date(2005, 1, 1), 4383]}

# Note on CAX - additional dataset from met office only has one location, not CAX06 and CAX04

# start date and number of daily time steps for lba obs files - CARBON
site_dates_lba = {
"BAN": [dt.date(2004, 1, 1), 1096],
"FNS": [dt.date(1999, 1, 1), 1096],
"K34": [dt.date(2003, 1, 1), 1461],
"K67": [dt.date(2001, 12, 31), 1462],
"K77": [dt.date(2001, 1, 1), 1826],
"K83": [dt.date(2001, 1, 1), 1095],
"PDG": [dt.date(2002, 1, 1), 730],
"RJA": [dt.date(2000, 1, 1), 1096],
"CAX": [dt.date(1999, 1, 1), 731]}

""" Not present in LBA OBS:
"CAX04": [dt.date(2005, 1, 1), 4383],
"CAX06": [dt.date(2005, 1, 1), 4383],
"KEN01": [dt.date(2005, 1, 1), 4383],
"KEN02": [dt.date(2005, 1, 1), 4383],
"TAM05": [dt.date(2005, 1, 1), 4383],
"TAM06": [dt.date(2005, 1, 1), 4383],
"NVX": [dt.date(2005, 1, 1), 4383]}
"""

if obs_data=="lba":
    path2csv="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/lba_obs/"
    path2nc="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/lba_obs/nc/"
elif obs_data=="gem":
    path2csv=="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/GEM_Sophie/"
    # special case here:


# inland output variables
out_vars = ["swnet", "lwnet", "qle", "qh", "qg"] # EDIT THIS, add variable names, name them as you like

lba_obs_out_vars = ["NEE","NEEf","NEE_model","Re_5day_ust_Sco2_LUT","gpp_gb","par_fill","VPD","mrs"]
lba_obs_in_vars = ["NEE","NEEf","NEE_model","Re_5day_ust_Sco2_LUT","GEP_model","par_fill","VPD","mrs"]

##########################################################################################################

def generate_dates_from_lba_obs_csv(dataframe_lba_obs_file):
    """
    TODO: (CURRENTLY DONE MANUALLY ABOVE)

    Generates a dictionary of sites:datetime.date objects from the LBA site csv file.

    Input:
        Pandas data frame with a single lba obs site loaded.

    Return:
        Dictionary of the sites:datetime.date objects
    """
    df = dataframe_lba_obs_file
    # Each dt.date object needs a start year, month and day
    # TODOsite_code = df.
    #TODOstart_year
    
def extract_carbon_lba_obs_csv(filename, site):
    """
    Takes a filename of a csv file for lba obs data and returns the file as a pandas dataframe.
    """
    print "##############################################"
    print "Extracting LBS OBS data from csv file at " + site

    CARBON_HEADERS = ["Year_LBAMIP","DoY_LBAMIP","Hour_LBAMIP","NEE","NEEf","NEE_model","Re_5day_ust_Sco2_LUT","GEP_model","par_fill","VPD","mrs"]

    import pandas as pd
    carbon_lba_obs_df = pd.read_csv(filename, delimiter=',', header=None, names=CARBON_HEADERS)
    return carbon_lba_obs_df

def write_lba_obs_netcdf(obs_output, site, outfilename):
    """
    Takes a pandas dataframe with the obs data a writes it to netcdf, for a given site
    """

    print "Writing LBA OBS output to netcdf with relevant metadata/attributes"

    data = Dataset(outfilename, 'w', format='NETCDF4')

    lats = [sitelocs[site][0]]
    lons = [sitelocs[site][1]]

    latdim=data.createDimension('lat', len(lats))
    londim=data.createDimension('lon', len(lons))
    tdim=data.createDimension('time', None) # record, or unlimited dimension
	
	# variables
    latitudes=data.createVariable('lat', 'f4', ('lat',))
    longitudes=data.createVariable('lon', 'f4', ('lon',))
    time=data.createVariable('time', np.float64, ('time',))

    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'

    latitudes.title = 'Latitude'
    longitudes.title = 'Longitude'

    latitudes.actual_min = min(lats)
    latitudes.actual_max = max(lats)

    longitudes.actual_min = min(lons)
    longitudes.actual_max = max(lons)

    latitudes[:] = lats
    longitudes[:] = lons

    print "Creating metadata"

    # add values to time variable

    # For lba this is daily! (Not hourly)
    #times = [datetime.datetime.combine(site_dates[site][0], datetime.time()) + datetime.timedelta(hours=hour) for hour in range(site_dates[site][1]*24)]
 
    #times = [datetime.date(site_dates[site][0] + datetime.timedelta(days=day) for day in range(site_dates[site][1])]

    # For daily timeseries, you can do this:
    times = [datetime.datetime.combine(site_dates_lba[site][0], datetime.time()) + datetime.timedelta(days=day) for day in range(site_dates_lba[site][1])]
    #times = [dt.datetime.combine(dates[i], dt.time()) for i, date in enumerate(dates)]

    time.units = 'days since 1850-01-01 00:00:00.0'  # CHECK THIS MATCHES WHAT IS IN YOUR INPUT DATA!
    time.calendar = 'gregorian'
    time[:] = date2num(times, time.units, calendar=time.calendar)

    print "Writing data to file for "+str(len(lats))+" point(s)"

    for in_var in lba_obs_in_vars:
        print "Creating netcdf variable: " + in_var
        dataout=data.createVariable(in_var, 'f4', ('time', 'lat', 'lon',),fill_value=-9999.)

		# EDIT THIS, add metadata for other variables here
        if in_var=="GEP_model":  # GPP_GB alias	
            dataout.units = "g.m-2.d-1"
            dataout.long_name = "Gridbox GPP"
            dataout.title = "GPP_GB"
        if in_var=="NEE_model": 
            dataout.units = "g.m-2.d-1"
            dataout.long_name = "Gridbox NEE"
        else:
            print "No data variables to write...? \n Check your netcdf file"

        print "Writing "+ in_var +" data to file..."
        print obs_output[in_var]
        dataout[:,0,0] = obs_output[in_var].values

		
        # min/max values
        dataout.actual_min = np.min(dataout[:,0,0])
        dataout.actual_max = np.max(dataout[:,0,0])

    data.close()
    print "***SUCCESS writing to file***"	

##########################################################################################################
# extract data from all sites and write to netcdf
for site in sites:

    print "PROCESSING SITE: " + site
    # Extract the data from the site csv
	#inland_outputs = extract_csv_data(path2csv+"/"+i+"-single_point-output.csv",i)
    lba_obs_outputs = extract_carbon_lba_obs_csv(path2csv+"/" + site + "day-carbon.dat", site)

    print "LBA OBS dataframe summary:"
    print lba_obs_outputs
	
    # Concatenate an output file csv
    outfile = path2nc+"/"+site+"-carbon_LBA_site_point.nc"

    # Write the variable to netcdf format
    write_lba_obs_netcdf(lba_obs_outputs, site, outfile)

##########################################################################################################


















