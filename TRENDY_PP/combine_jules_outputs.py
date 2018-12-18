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
varlist=["gpp", "npp", "cveg"]

startyear=2000
startmonth=1

endyear=2017
endmonth=12

nmonths=216

path2data="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/TRENDY/MODEL_DATA/JULES-ES.1p0.vn5.1.500.CRUJRA1.365.S3.Monthly.2D"

##########################################################################################################
def extract_data(modelvar):
	
	if len(str(startmonth))==1:
		beginmnth="0"+str(startmonth)
	elif len(str(startmonth))==2:
		beginmnth=str(startmonth)

	if len(str(endmonth))==1:
		endmnth="0"+str(endmonth)
	elif len(str(endmonth))==2:
		endmnth=str(endmonth)
	
	julesfile=path2data+"/"+modelvar+"/JULES-ES.1p0.vn5.1.500.CRUJRA1.365.S3.Monthly.2D."+str(startyear)+beginmnth+"-"+str(endyear)+endmnth+"."+modelvar+".nc"

	print "Extracting JULES output from "+julesfile

	data = Dataset(julesfile, "r")


	if modelvar=="gpp" or modelvar=="npp" or modelvar=="cveg":
		var_data = np.array(data.variables[modelvar][:])

	

	lats = np.array(data.variables["latitude"][:])
	lons = np.array(data.variables["longitude"][:])
	
	data.close()
	
	return var_data, lats, lons 


def write2netCDF(moutput, lats, lons, outfilename):

	print "Writing JULES output (with metadata, evaluation tool ready) to file"
	print "Creating file "+outfilename	

	data = Dataset(outfilename, 'w', format='NETCDF4')
	
	# dimensions
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

	print "Writing data to file for "+str(len(lats))+" point(s)"

	# add values to time variable
	times = [datetime.datetime.combine(dt.date(2000,1,1), datetime.time()) + datetime.timedelta(hours=month) for month in range(nmonths)]
	time.units = 'days since 1850-01-01 00:00:00.0'
	time.calendar = 'gregorian'
	time[:] = date2num(times, time.units, calendar=time.calendar)	

	for i in varlist:

		dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=9.96920996839e+36)

		# EDIT THIS, add metadata for other variables here		
		if i=="gpp":	
			dataout.units = "g.m-2.d-1"
			dataout.long_name = "Gridbox GPP"
		elif i=="npp":		
			dataout.units = "g.m-2.d-1"
			dataout.long_name = "Gridbox net primary productivity prior to N limitation"
		elif i=="cveg":
			dataout.units = "g.m-2"			
			dataout.long_name = "Gridbox mean vegetation carbon at end of model timestep"

		print "Writing "+i+" data to file..."
		
		dataout[:,:,:] = moutput[i][:]

		if i=="gpp" or i=="npp":
			dataout[:,:,:]*=1000*86400 # kg.m-2.s-1 -> g.m-2.day-1
		elif i=="cveg":
			dataout[:,:,:]*=1000 # kg.m-2 -> g.m-2		

		# min/max values
		dataout.actual_min = np.min(dataout[:,:,:])
		dataout.actual_max = np.max(dataout[:,:,:])

	data.close()
	print "***SUCCESS writing to file***"	


##########################################################################################################
modeloutput={}

for outputvar in varlist:
	print "Extracting "+outputvar+" data"	
	
	model_tmp, lat_tmp, lon_tmp=extract_data(outputvar)
	modeloutput[outputvar]=model_tmp

outfile = path2data+"/JULES-ES.1p0.vn5.1.500.CRUJRA1.365.S3.Monthly.2D.ILAMB.nc"
write2netCDF(modeloutput, lat_tmp, lon_tmp, outfile)

	


















































































