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
# Met data can be local (LBA, GEM) or global (CRUJRA)
met_data="global"

# add sites to this array for model runs
#sites=["BAN", "FNS", "K34", "K67", "K77", "K83", "PDG", "RJA", "CAX04", "CAX06", "KEN01", "KEN02", "TAM05", "TAM06", "NVX"]
if met_data=="local":
    sites=["BAN","BAN-veg2","BAN-veg9","CAX04","CAX06","FNS","K34","K67","K77","K83","NVX","PDG","RJA","TAM05","TAM06"] # EDIT THIS, add site names to this array
elif met_data=="global":
    sites=["BAN","CAX04","CAX06","FNS","K34","K67","K77","K83","KEN01","NVX","PDG","RJA","TAM05","TAM06"] # EDIT THIS, add site names to this array

# site: [lat, lon]
sitelocs = {"CAX04": [-1.716, -51.457], "CAX06": [-1.737, -51.462], "KEN01": [-16.016, -62.730], "KEN02": [-16.016, -62.730],
"TAM05": [-12.831, -69.271], "TAM06": [-12.839, -69.296], "NVX": [-14.70, -52.35], "BAN": [-9.82, -50.15], "FNS": [-10.76, -62.36],
"K34": [-2.61, -60.21], "K67": [-2.85, -54.97], "K77": [-3.01, -54.54], "K83": [-3.05, -54.93],
"PDG": [-21.62, -47.63], "RJA": [-10.08, -62.36], "BAN-veg2": [-9.82, -50.15], "BAN-veg9": [-9.82, -50.15]}

# start date and number of daily time steps
site_dates = {"BAN": [dt.date(2004, 1, 1), 1095], "FNS": [dt.date(1999, 1, 1), 1096], "K34": [dt.date(2003, 1, 1), 1096],\
 "K67": [dt.date(2002, 1, 1), 1096], "K77": [dt.date(2001, 1, 1), 1826], "K83": [dt.date(2001, 1, 1), 1095],\
 "PDG": [dt.date(2002, 1, 1), 730], "RJA": [dt.date(2000, 1, 1), 1096], "CAX04": [dt.date(2005, 1, 1), 4383],\
 "CAX06": [dt.date(2005, 1, 1), 4383], "KEN01": [dt.date(2005, 1, 1), 4383], "KEN02": [dt.date(2005, 1, 1), 4383],\
 "TAM05": [dt.date(2005, 1, 1), 4383], "TAM06": [dt.date(2005, 1, 1), 4383], "NVX": [dt.date(2005, 1, 1), 4383],\
 "BAN-veg2": [dt.date(2004, 1, 2), 1095], "BAN-veg9": [dt.date(2004, 1, 2), 1095]}

#site_dates = {"BAN": [dt.date(2004, 1, 2), 1033], "FNS": [dt.date(1999, 1, 2), 1094], "K34": [dt.date(2003, 1, 2), 1017], "K67": [dt.date(2002, 1, 2), 685],\
#"K77": [dt.date(2001, 1, 2), 1824], "K83": [dt.date(2001, 1, 2), 952], "PDG": [dt.date(2002, 1, 2), 728], "RJA": [dt.date(2000, 2, 3), 953],\
#"CAX04": [dt.date(2005, 1, 1), 4383], "CAX06": [dt.date(2005, 1, 1), 4383], "KEN01": [dt.date(2005, 1, 1), 4383], "KEN02": [dt.date(2005, 1, 1), 4383],\
#"TAM05": [dt.date(2005, 1, 1), 4383], "TAM06": [dt.date(2005, 1, 1), 4383], "NVX": [dt.date(2005, 1, 1), 4383], "BAN-veg2": [dt.date(2004, 1, 2), 1033],\
#"BAN-veg9": [dt.date(2004, 1, 2), 1033]}

if met_data=="local":
	path2csv="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/model_outputs/INLAND/site_runs/met/local/csv"
	path2nc="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/model_outputs/INLAND/site_runs/met/local/nc"
#	path2csv="/Users/aline/Modelos/INLAND-OFFLINE-MRSO-GRIDxSINGLE/inland-offline-single-withCvR/output/LOCAL/csv"
#	path2nc="/Users/aline/Modelos/INLAND-OFFLINE-MRSO-GRIDxSINGLE/inland-offline-single-withCvR/output/LOCAL/nc"
elif met_data=="global":
	path2csv="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/model_outputs/INLAND/site_runs/met/global/csv"
	path2nc="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/model_outputs/INLAND/site_runs/met/global/nc"
#	path2csv="/Users/aline/Modelos/INLAND-OFFLINE-MRSO-GRIDxSINGLE/inland-offline-single-withCvR/output/CRUJRA/csv"
#	path2nc="/Users/aline/Modelos/INLAND-OFFLINE-MRSO-GRIDxSINGLE/inland-offline-single-withCvR/output/CRUJRA/nc"

# inland output variables
#out_vars = ["swnet", "lwnet", "qle", "qh", "qg"] #  EDIT THIS, add variable names, name them as you like
out_vars = ["Swnet","Lwnet","Qle","Qh","Qg","DelCanHeat","Evap","Qs","Qsb","DelSoilMoist","VegT","BaresoilT","Albedo","fPAR",\
            "SoilMoist1","SoilMoist2","SoilMoist3","SoilMoist4","SoilMoist5","SoilMoist6","SoilTemp1","SoilTemp2","SoilTemp3",\
            "SoilTemp4","SoilTemp5","SoilTemp6","SoilWet","Ecanop","Tveg","Esoil","CanopInt","GPP","NPP","NEE","AutoResp","HeteroResp",\
            "AGBNPP", "LeafLitterFall","RootTurnover", "WoodTurnover","cbioltot","cbiortot","cbiowtot","clitll","clitlm","clitls",\
            "clitrl", "clitrm","clitrs","clitwl","clitwm","clitws","AGB","TotLivBiom","AGB_Wood","LAI",\
            "Tair","Qair","Wind","Rainf","Psurf","SWdown","LWdown","CO2air"] #  EDIT THIS, add variable names, name them as you like

##########################################################################################################
def extract_csv_data(filename,site):

	print "#############################################################################################"
	print "Extracting INLAND data from csv file at "+site

	# store data in a dictionary
	inland_data = {}

	data = np.loadtxt(filename, skiprows=1, dtype=np.str)
	
	###########################################################################################
	# EDIT THIS
	# The column number (starting at 0) is used to extract data from the data array. 	
	###########################################################################################
	if met_data=="local":
		inland_data["Swnet"] = np.array([float(x) for x in data[:,3]])
		inland_data["Lwnet"] = np.array([float(x) for x in data[:,4]])
		inland_data["Qle"] = np.array([float(x) for x in data[:,5]])
		inland_data["Qh"] = np.array([float(x) for x in data[:,6]])
		inland_data["Qg"] = np.array([float(x) for x in data[:,7]])
		inland_data["DelCanHeat"] = np.array([float(x) for x in data[:,8]])
		inland_data["Evap"] = np.array([float(x) for x in data[:,10]])
		inland_data["Qs"] = np.array([float(x) for x in data[:,11]])
		inland_data["Qsb"] = np.array([float(x) for x in data[:,13]])
		inland_data["DelSoilMoist"] = np.array([float(x) for x in data[:,15]])
		inland_data["VegT"] = np.array([float(x) for x in data[:,17]])
		inland_data["BaresoilT"] = np.array([float(x) for x in data[:,18]])
		inland_data["Albedo"] = np.array([float(x) for x in data[:,20]])
		inland_data["fPAR"] = np.array([float(x) for x in data[:,22]])
		inland_data["SoilMoist1"] = np.array([float(x) for x in data[:,23]])
		inland_data["SoilMoist2"] = np.array([float(x) for x in data[:,24]])
		inland_data["SoilMoist3"] = np.array([float(x) for x in data[:,25]])
		inland_data["SoilMoist4"] = np.array([float(x) for x in data[:,26]])
		inland_data["SoilMoist5"] = np.array([float(x) for x in data[:,27]])
		inland_data["SoilMoist6"] = np.array([float(x) for x in data[:,28]])
		inland_data["SoilTemp1"] = np.array([float(x) for x in data[:,29]])
		inland_data["SoilTemp2"] = np.array([float(x) for x in data[:,30]])
		inland_data["SoilTemp3"] = np.array([float(x) for x in data[:,31]])
		inland_data["SoilTemp4"] = np.array([float(x) for x in data[:,32]])
		inland_data["SoilTemp5"] = np.array([float(x) for x in data[:,33]])
		inland_data["SoilTemp6"] = np.array([float(x) for x in data[:,34]])
		inland_data["SoilWet"] = np.array([float(x) for x in data[:,35]])
		inland_data["Ecanop"] = np.array([float(x) for x in data[:,36]])
		inland_data["Tveg"] = np.array([float(x) for x in data[:,37]])
		inland_data["Esoil"] = np.array([float(x) for x in data[:,38]])
		inland_data["RootMoist"] = np.array([float(x) for x in data[:,40]])
		inland_data["CanopInt"] = np.array([float(x) for x in data[:,41]])
		inland_data["GPP"] = np.array([float(x)*1000*86400 for x in data[:,42]]) # kgC.m-2.s-1 --> gC.m-2.day-1
		inland_data["NPP"] = np.array([float(x)*1000*86400 for x in data[:,43]]) # kgC.m-2.s-1 --> gC.m-2.day-1
		inland_data["NEE"] = np.array([float(x)*1000*86400 for x in data[:,44]]) # kgC.m-2.s-1 --> gC.m-2.day-1
		inland_data["AutoResp"] = np.array([float(x)*1000*86400 for x in data[:,45]]) # kgC.m-2.s-1 --> gC.m-2.day-1
		inland_data["HeteroResp"] = np.array([float(x)*1000*86400 for x in data[:,46]]) # kgC.m-2.s-1 --> gC.m-2.day-1
		inland_data["AGBNPP"] = np.array([float(x)*1000 for x in data[:,47]]) # kgC.m-2.year-1 --> gC.m-2.year-1
		inland_data["Rootlitter"] = np.array([float(x) for x in data[:,48]])
		inland_data["LeafLitterFall"] = np.array([float(x)*1000 for x in data[:,49]]) # kgC.m-2.year-1 --> gC.m-2.year-1
		inland_data["RootTurnover"] = np.array([float(x)*1000 for x in data[:,50]]) # kgC.m-2.year-1 --> gC.m-2.year-1
		inland_data["WoodTurnover"] = np.array([float(x)*1000 for x in data[:,51]]) # kgC.m-2.year-1 --> gC.m-2.year-1
		inland_data["cbioltot"] = np.array([float(x)*1000 for x in data[:,52]]) # kgC.m-2 --> gC.m-2
		inland_data["cbiortot"] = np.array([float(x)*1000 for x in data[:,53]]) # kgC.m-2 --> gC.m-2
		inland_data["cbiowtot"] = np.array([float(x)*1000 for x in data[:,54]]) # kgC.m-2 --> gC.m-2
		inland_data["clitll"] = np.array([float(x)*1000 for x in data[:,55]]) # kgC.m-2 --> gC.m-2
		inland_data["clitlm"] = np.array([float(x)*1000 for x in data[:,56]]) # kgC.m-2 --> gC.m-2
		inland_data["clitls"] = np.array([float(x)*1000 for x in data[:,57]]) # kgC.m-2 --> gC.m-2
		inland_data["clitrl"] = np.array([float(x)*1000 for x in data[:,58]]) # kgC.m-2 --> gC.m-2
		inland_data["clitrm"] = np.array([float(x)*1000 for x in data[:,59]]) # kgC.m-2 --> gC.m-2
		inland_data["clitrs"] = np.array([float(x)*1000 for x in data[:,60]]) # kgC.m-2 --> gC.m-2
		inland_data["clitwl"] = np.array([float(x)*1000 for x in data[:,61]]) # kgC.m-2 --> gC.m-2
		inland_data["clitwm"] = np.array([float(x)*1000 for x in data[:,62]]) # kgC.m-2 --> gC.m-2
		inland_data["clitws"] = np.array([float(x)*1000 for x in data[:,63]]) # kgC.m-2 --> gC.m-2
		inland_data["AGB"] = np.array([float(x)*1000 for x in data[:,64]]) # kgC.m-2 --> gC.m-2
		inland_data["TotLivBiom"] = np.array([float(x)*1000 for x in data[:,65]]) # kgC.m-2 --> gC.m-2
		inland_data["AGB_Wood"] = np.array([float(x)*1000 for x in data[:,66]]) # kgC.m-2 --> gC.m-2
		inland_data["LAI"] = np.array([float(x) for x in data[:,67]])
		inland_data["Tair"] = np.array([float(x) for x in data[:,68]])
		inland_data["Qair"] = np.array([float(x) for x in data[:,69]])
		inland_data["Wind"] = np.array([float(x) for x in data[:,70]])
		inland_data["Rainf"] = np.array([float(x) for x in data[:,71]])
		inland_data["Psurf"] = np.array([float(x) for x in data[:,72]])
		inland_data["SWdown"] = np.array([float(x) for x in data[:,73]])
		inland_data["LWdown"] = np.array([float(x) for x in data[:,74]])
		inland_data["CO2air"] = np.array([float(x) for x in data[:,75]])

	elif met_data=="global":
		inland_data["Swnet"] = np.array([float(x) for x in data[:,3]])
		inland_data["Lwnet"] = np.array([float(x) for x in data[:,4]])
		inland_data["Qle"] = np.array([float(x) for x in data[:,5]])
		inland_data["Qh"] = np.array([float(x) for x in data[:,6]])
		inland_data["Qg"] = np.array([float(x) for x in data[:,7]])
		inland_data["DelCanHeat"] = np.array([float(x) for x in data[:,8]])
		inland_data["Evap"] = np.array([float(x) for x in data[:,9]])
		inland_data["Qs"] = np.array([float(x) for x in data[:,10]])
		inland_data["Qsb"] = np.array([float(x) for x in data[:,11]])
		inland_data["DelSoilMoist"] = np.array([float(x) for x in data[:,12]])
		inland_data["VegT"] = np.array([float(x) for x in data[:,13]])
		inland_data["BaresoilT"] = np.array([float(x) for x in data[:,14]])
		inland_data["Albedo"] = np.array([float(x) for x in data[:,15]])
		inland_data["fPAR"] = np.array([float(x) for x in data[:,16]])
		inland_data["SoilMoist1"] = np.array([float(x) for x in data[:,17]])
		inland_data["SoilMoist2"] = np.array([float(x) for x in data[:,18]])
		inland_data["SoilMoist3"] = np.array([float(x) for x in data[:,19]])
		inland_data["SoilMoist4"] = np.array([float(x) for x in data[:,20]])
		inland_data["SoilMoist5"] = np.array([float(x) for x in data[:,21]])
		inland_data["SoilMoist6"] = np.array([float(x) for x in data[:,22]])
		inland_data["SoilTemp1"] = np.array([float(x) for x in data[:,23]])
		inland_data["SoilTemp2"] = np.array([float(x) for x in data[:,24]])
		inland_data["SoilTemp3"] = np.array([float(x) for x in data[:,25]])
		inland_data["SoilTemp4"] = np.array([float(x) for x in data[:,26]])
		inland_data["SoilTemp5"] = np.array([float(x) for x in data[:,27]])
		inland_data["SoilTemp6"] = np.array([float(x) for x in data[:,28]])
		inland_data["SoilWet"] = np.array([float(x) for x in data[:,29]])
		inland_data["Ecanop"] = np.array([float(x) for x in data[:,30]])
		inland_data["Tveg"] = np.array([float(x) for x in data[:,31]])
		inland_data["Esoil"] = np.array([float(x) for x in data[:,32]])
		inland_data["CanopInt"] = np.array([float(x) for x in data[:,33]])
		inland_data["GPP"] = np.array([float(x)*1000*86400 for x in data[:,33]]) # kgC.m-2.s-1 --> gC.m-2.day-1
		inland_data["NPP"] = np.array([float(x)*1000*86400 for x in data[:,34]]) # kgC.m-2.s-1 --> gC.m-2.day-1
		inland_data["NEE"] = np.array([float(x)*1000*86400 for x in data[:,35]]) # kgC.m-2.s-1 --> gC.m-2.day-1
		inland_data["AutoResp"] = np.array([float(x)*1000*86400 for x in data[:,36]]) # kgC.m-2.s-1 --> gC.m-2.day-1
		inland_data["HeteroResp"] = np.array([float(x)*1000*86400 for x in data[:,37]]) # kgC.m-2.s-1 --> gC.m-2.day-1
		inland_data["AGBNPP"] = np.array([float(x)*1000 for x in data[:,39]]) # kgC.m-2.year-1 --> gC.m-2.year-1
		inland_data["LeafLitterFall"] = np.array([float(x)*1000 for x in data[:,40]]) # kgC.m-2.year-1 --> gC.m-2.year-1
		inland_data["RootTurnover"] = np.array([float(x)*1000 for x in data[:,41]]) # kgC.m-2.year-1 --> gC.m-2.year-1
		inland_data["WoodTurnover"] = np.array([float(x)*1000 for x in data[:,42]]) # kgC.m-2.year-1 --> gC.m-2.year-1
		inland_data["cbioltot"] = np.array([float(x)*1000 for x in data[:,43]]) # kgC.m-2 --> gC.m-2
		inland_data["cbiortot"] = np.array([float(x)*1000 for x in data[:,44]]) # kgC.m-2 --> gC.m-2
		inland_data["cbiowtot"] = np.array([float(x)*1000 for x in data[:,45]]) # kgC.m-2 --> gC.m-2
		inland_data["clitll"] = np.array([float(x)*1000 for x in data[:,46]]) # kgC.m-2 --> gC.m-2
		inland_data["clitlm"] = np.array([float(x)*1000 for x in data[:,47]]) # kgC.m-2 --> gC.m-2
		inland_data["clitls"] = np.array([float(x)*1000 for x in data[:,48]]) # kgC.m-2 --> gC.m-2
		inland_data["clitrl"] = np.array([float(x)*1000 for x in data[:,49]]) # kgC.m-2 --> gC.m-2
		inland_data["clitrm"] = np.array([float(x)*1000 for x in data[:,50]]) # kgC.m-2 --> gC.m-2
		inland_data["clitrs"] = np.array([float(x)*1000 for x in data[:,51]]) # kgC.m-2 --> gC.m-2
		inland_data["clitwl"] = np.array([float(x)*1000 for x in data[:,52]]) # kgC.m-2 --> gC.m-2
		inland_data["clitwm"] = np.array([float(x)*1000 for x in data[:,53]]) # kgC.m-2 --> gC.m-2
		inland_data["clitws"] = np.array([float(x)*1000 for x in data[:,54]]) # kgC.m-2 --> gC.m-2
		inland_data["AGB"] = np.array([float(x)*1000 for x in data[:,55]]) # kgC.m-2 --> gC.m-2
		inland_data["TotLivBiom"] = np.array([float(x)*1000 for x in data[:,56]]) # kgC.m-2 --> gC.m-2
		inland_data["AGB_Wood"] = np.array([float(x)*1000 for x in data[:,57]]) # kgC.m-2 --> gC.m-2
		inland_data["LAI"] = np.array([float(x) for x in data[:,58]])
		inland_data["Tair"] = np.array([float(x) for x in data[:,59]])
		inland_data["Qair"] = np.array([float(x) for x in data[:,60]])
		inland_data["Wind"] = np.array([float(x) for x in data[:,61]])
		inland_data["Rainf"] = np.array([float(x) for x in data[:,62]])
		inland_data["Psurf"] = np.array([float(x) for x in data[:,63]])
		inland_data["SWdown"] = np.array([float(x) for x in data[:,64]])
		inland_data["LWdown"] = np.array([float(x) for x in data[:,65]])
		inland_data["CO2air"] = np.array([float(x) for x in data[:,66]])


	###########################################################################################

	return inland_data


def write2netCDF(moutput, site, outfilename):

	print "Writing INLAND output (with metadata, evaluation tool, ILAMB, ready) to file"
	
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
	if met_data=="local":
		times = [datetime.datetime.combine(site_dates[site][0], datetime.time()) + datetime.timedelta(hours=hour) for hour in range(site_dates[site][1]*24)]
	elif met_data=="global":
	 	times = [datetime.datetime.combine(dt.date(1901, 1, 1), datetime.time()) + datetime.timedelta(hours=hour) for hour in range(42734*24)]

    	time.units = 'hours since 1850-01-01 00:00:00.0'
    	time.calendar = 'gregorian'
	time[:] = date2num(times, time.units, calendar=time.calendar)

	print "Writing data to file for "+str(len(lats))+" point(s)"

	for i in out_vars:

		#dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)

		# EDIT THIS, add metadata for other variables here
		if i=="Swnet":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "W.m-2"
			dataout.long_name = "Net shortwave radiation"
			dataout.title = "Swnet"
		elif i=="Lwnet":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "W.m-2"
			dataout.long_name = "Net longwave radiation"
			dataout.title = "Lwnet"
		elif i=="Qle":	
			dataout=data.createVariable("latent_heat", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "W.m-2"
			dataout.long_name = "Latent heat"
			dataout.title = "Qle"
		elif i=="Qh":
			dataout=data.createVariable("ftl_gb", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "W.m-2"
			dataout.long_name = "Sensible heat"
			dataout.title = "Qh"
		elif i=="Qg":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "W.m-2"
			dataout.long_name = "Soil heat flux"
			dataout.title = "Qg"
		elif i=="DelCanHeat":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "W.m-2"
			dataout.long_name = "Change incanopy heat storage"
			dataout.title = "DelCanHeat"
		elif i=="Evap":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "mm.s-1"
			dataout.long_name = "Evapotranspiration"
			dataout.title = "Evap"
		elif i=="Qs":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "mm.s-1"
			dataout.long_name = "Surface runoff"
			dataout.title = "Qs"
		elif i=="Qsb":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "mm.s-1"
			dataout.long_name = "Subsurface runoff"
			dataout.title = "Qsb"
		elif i=="DelSoilMoist":
			dataout=data.createVariable("smc_tot", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg.m-2"
			dataout.long_name = "Change in soil moisture"
			dataout.title = "DelSoilMoist"
		elif i=="VegT":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "K"
			dataout.long_name = "Vegetation canopy temperature"
			dataout.title = "VegT"
		elif i=="BaresoilT":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "K"
			dataout.long_name = "Soil temperature in the first soil layer"
			dataout.title = "BaresoilT"
		elif i=="Albedo":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "fraction"
			dataout.long_name = "Surface albedo"
			dataout.title = "Albedo"
		elif i=="fPAR":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "fraction"
			dataout.long_name = "Asorbed fraction of PAR"
			dataout.title = "fPAR"
		elif i=="SoilMoist1":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg.m-2"
			dataout.long_name = "Soils moisture - layer 1"
			dataout.title = "SoilMoist1"
		elif i=="SoilMoist2":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg.m-2"
			dataout.long_name = "Soils moisture - layer 2"
			dataout.title = "SoilMoist2"
		elif i=="SoilMoist3":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg.m-2"
			dataout.long_name = "Soils moisture - layer 3"
			dataout.title = "SoilMoist3"
		elif i=="SoilMoist4":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg.m-2"
			dataout.long_name = "Soils moisture - layer 4"
			dataout.title = "SoilMoist4"
		elif i=="SoilMoist5":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg.m-2"
			dataout.long_name = "Soils moisture - layer 5"
			dataout.title = "SoilMoist5"
		elif i=="SoilMoist6":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg.m-2"
			dataout.long_name = "Soils moisture - layer 6"
			dataout.title = "SoilMoist6"
		elif i=="SoilTemp1":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "K"
			dataout.long_name = "Soils temperature - layer 1"
			dataout.title = "SoilTemp1"
		elif i=="SoilTemp2":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "K"
			dataout.long_name = "Soils temperature - layer 2"
			dataout.title = "SoilTemp2"
		elif i=="SoilTemp3":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "K"
			dataout.long_name = "Soils temperature - layer 3"
			dataout.title = "SoilTemp3"
		elif i=="SoilTemp4":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "K"
			dataout.long_name = "Soils temperature - layer 4"
			dataout.title = "SoilTemp4"
		elif i=="SoilTemp5":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "K"
			dataout.long_name = "Soils temperature - layer 5"
			dataout.title = "SoilTemp5"
		elif i=="SoilTemp6":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "K"
			dataout.long_name = "Soils temperature - layer 6"
			dataout.title = "SoilTemp6"
		elif i=="SoilWet":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "fraction?"
			dataout.long_name = "Total soil wetness"
			dataout.title = "SoilWet"
		elif i=="Ecanop":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "mm.s-1"
			dataout.long_name = "Total evaporation rate from all intercepted h2o"
			dataout.title = "Ecanop"
		elif i=="Tveg":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "mm.s-1"
			dataout.long_name = "Total transpiration rate from all vegetation canopies"
			dataout.title = "Tveg"
		elif i=="Esoil":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "mm.s-1"
			dataout.long_name = "Total evaporation rate from surface (snow/soil)"
			dataout.title = "Esoil"
		elif i=="CanopInt":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "mm?"
			dataout.long_name = "Total canopy water storage"
			dataout.title = "CanopInt"
		elif i=="GPP":
			dataout=data.createVariable("gpp_gb", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)
			dataout.units = "0.012.(mol_CO2.m-2.s-1)"
			dataout.long_name = "Gross primary productivity"
			dataout.title = "GPP"
		elif i=="NPP":
			dataout=data.createVariable("npp_gb", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "0.012.(mol_CO2.m-2.s-1)"
			dataout.long_name = "Net primary productivity"
			dataout.title = "NPP"
		elif i=="NEE":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "0.012.(mol_CO2.m-2.s-1)"
			dataout.long_name = "Net ecosystem exchange"
			dataout.title = "NEE"
		elif i=="AutoResp":
			dataout=data.createVariable("resp_p_gb", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "0.012.(mol_CO2.m-2.s-1)"
			dataout.long_name = "Autotrophic respiration"
			dataout.title = "AutoResp"
		elif i=="HeteroResp":
			dataout=data.createVariable("resp_s_gb", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "0.012.(mol_CO2.m-2.s-1)"
			dataout.long_name = "Heterotrophic respiration"
			dataout.title = "HeteroResp"
		elif i=="AGBNPP":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2.year-1"
			dataout.long_name = "Annual above-ground npp for ecosystem"
			dataout.title = "AGBNPP"
		elif i=="LeafLitterFall":
			dataout=data.createVariable("leaf_litC", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2.year-1"
			dataout.long_name = "Annual leaf litter fall"
			dataout.title = "LeafLitterFall"
		elif i=="RootTurnover":
			dataout=data.createVariable("root_litC", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2.year-1"
			dataout.long_name = "Annual fine root turnover"
			dataout.title = "RootTurnover"
		elif i=="WoodTurnover":
			dataout=data.createVariable("wood_litC", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2.year-1"
			dataout.long_name = "Annual woody turnover"
			dataout.title = "WoodTurnover"
		elif i=="cbioltot":
			dataout=data.createVariable("leafC", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Sum in PFTs carbon in leaf biomass pool"
			dataout.title = "cbioltot"
		elif i=="cbiortot":
			dataout=data.createVariable("rootC", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Sum in PFTs carbon in fine root biomass pool"
			dataout.title = "cbiortot"
		elif i=="cbiowtot":
			dataout=data.createVariable("woodC", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Sum in PFTs carbon in woody biomass pool"
			dataout.title = "cbiowtot"
		elif i=="clitll":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Carbon in leaf litter pool - lignin"
			dataout.title = "clitll"
		elif i=="clitlm":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Carbon in leaf litter pool - metabolic"
			dataout.title = "clitlm"
		elif i=="clitls":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Carbon in leaf litter pool - structural"
			dataout.title = "clitls"
		elif i=="clitrl":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Carbon in fine root litter pool - lignin"
			dataout.title = "clitrl"
		elif i=="clitrm":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Carbon in fine root litter pool - metabolic"
			dataout.title = "clitrm"
		elif i=="clitrs":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Carbon in fine root litter pool - structural"
			dataout.title = "clitrs"
		elif i=="clitwl":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Carbon in woody litter pool - lignin"
			dataout.title = "clitwl"
		elif i=="clitwm":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Carbon in woody litter pool - metabolic"
			dataout.title = "clitwm"
		elif i=="clitws":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Carbon in woody litter pool - structural"
			dataout.title = "clitws"
		elif i=="AGB":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Above ground biomass"
			dataout.title = "AGB"
		elif i=="TotLivBiom":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Total living biomass"
			dataout.title = "TotLivBiom"
		elif i=="AGB_Wood":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg_C.m-2"
			dataout.long_name = "Above ground woody biomass"
			dataout.title = "AGB_Wood"
		elif i=="LAI":
			dataout=data.createVariable("lai", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = ""
			dataout.long_name = "Leaf area index"
			dataout.title = "LAI"
		elif i=="Tair":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "K"
			dataout.long_name = "Air temperature"
			dataout.title = "Tair"
		elif i=="Qair":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "kg.kg-1"
			dataout.long_name = "Specific humidity"
			dataout.title = "Qair"
		elif i=="Wind":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "m.s-1"
			dataout.long_name = "Wind"
			dataout.title = "Wind"
		elif i=="Rainf":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "mm.s-1"
			dataout.long_name = "Precipitation"
			dataout.title = "Rainf"
		elif i=="Psurf":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "Pa"
			dataout.long_name = "Surface pressure"
			dataout.title = "Psurf"
		elif i=="SWdown":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "W.m-2"
			dataout.long_name = "Shortwave incoming radiation"
			dataout.title = "SWdown"
		elif i=="LWdown":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "W.m-2"
			dataout.long_name = "Longwave incoming radiation"
			dataout.title = "LWdown"
		elif i=="CO2air":
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	
			dataout.units = "ppm?"
			dataout.long_name = "CO2 concentration"
			dataout.title = "CO2air"

		print "Writing "+i+" data to file..."

		#j=0
		#while j<len(time):
		dataout[:,0,0] = moutput[i][:]
		#	j+=1		

		# min/max values
		dataout.actual_min = min(dataout[:,0,0])
		dataout.actual_max = max(dataout[:,0,0])

	data.close()
	print "***SUCCESS writing to file***"	

##########################################################################################################
# extract data from all sites and write to netcdf
for i in sites:
	inland_outputs = extract_csv_data(path2csv+"/"+i+"-single_point-output.csv",i)
	
	outfile=path2nc+"/"+i+"-single_point-output_hourly.nc"

	write2netCDF(inland_outputs,i,outfile)

##########################################################################################################


















