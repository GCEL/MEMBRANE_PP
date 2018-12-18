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
from datetime import date
from dateutil.relativedelta import relativedelta

import pandas as pd
import xarray as xr

##########################################################################################################
# Can be lba, gem_lai, gem_npp (data from Doughtly et al 2015, Nature)
obs_data="gem_npp"

# temporal resolution of the output data
# LBA data is hourly.
# GEM LAI data is daily, monthly
# GEM NPP is monthly
timeres="monthly"

#sites=["BAN", "FNS", "K34", "K67", "K77", "K83", "PDG", "RJA", "CAX04", "CAX06", "KEN01", "KEN02", "TAM05", "TAM06", "NVX"]
# EDIT THIS, add site names to these arrays
if obs_data=="lba":
	sites=["BAN","FNS","CAX","K34","K67","K77","K83","PDG","RJA"]
elif obs_data=="gem_lai":
	sites=["KEN01", "KEN02", "TAM05", "TAM06", "CAX04", "CAX06"] # only have LAI for these GEM sites
elif obs_data=="gem_sm":
	sites=["KEN01"] # only have soil moisture for these GEM sites
elif obs_data=="gem_npp":
	sites=["CAX06"]





# site: [lat, lon]
sitelocs = {"CAX04": [-1.716, -51.457], "CAX06": [-1.737, -51.462], "KEN01": [-16.016, -62.730], "KEN02": [-16.016, -62.730],
"TAM05": [-12.831, -69.271], "TAM06": [-12.839, -69.296], "NVX": [-14.70, -52.35], "BAN": [-9.82, -50.15], "FNS": [-10.76, -62.36],
"K34": [-2.61, -60.21], "K67": [-2.85, -54.97], "K77": [-3.01, -54.54], "K83": [-3.05, -54.93],
"PDG": [-21.62, -47.63], "RJA": [-10.08, -62.36], "CAX": [-1.748, -51.454]}

# Note on CAX - lat lon taken from the site info from the LCB dataset:
# file:///exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/FLUX_LBA_ECO/CD32_BRAZIL_FLUX_NETWORK_1174/guide/CD32_Brazil_Flux_Network.html

# start date and number of daily time steps
site_dates = {"BAN": [dt.date(2004, 1, 2), 1033], "FNS": [dt.date(1999, 1, 2), 1094], "K34": [dt.date(2003, 1, 2), 1017], "K67": [dt.date(2002, 1, 2), 685],\
"K77": [dt.date(2001, 1, 2), 1824], "K83": [dt.date(2001, 1, 2), 952], "PDG": [dt.date(2002, 1, 2), 728], "RJA": [dt.date(2000, 2, 3), 953], "CAX": [dt.date(1999, 1 , 1), 731],\
"CAX04": [dt.date(2005, 1, 1), 4383], "CAX06": [dt.date(2005, 1, 1), 4383], "KEN01": [dt.date(2005, 1, 1), 4383], "KEN02": [dt.date(2005, 1, 1), 4383],\
"TAM05": [dt.date(2005, 1, 1), 4383], "TAM06": [dt.date(2005, 1, 1), 4383], "NVX": [dt.date(2005, 1, 1), 4383]}

# Note on CAX - additional dataset from met office only has one location, not CAX06 and CAX04

# start date and number of daily time steps for lba obs files - CARBON
site_dates = {
"BAN": [dt.date(2004, 1, 1), 1096],
"FNS": [dt.date(1999, 1, 1), 1096],
"K34": [dt.date(2003, 1, 1), 1461],
"K67": [dt.date(2001, 12, 31), 1462],
"K77": [dt.date(2001, 1, 1), 1826],
"K83": [dt.date(2001, 1, 1), 1095],
"PDG": [dt.date(2002, 1, 1), 730],
"RJA": [dt.date(2000, 1, 1), 1096],
"CAX": [dt.date(1999, 1, 1), 731],
"CAX04": [dt.date(2005, 1, 1), 4383],
"CAX06": [dt.date(2005, 1, 1), 4383],
"KEN01": [dt.date(2005, 1, 1), 4383],
"KEN02": [dt.date(2005, 1, 1), 4383],
"TAM05": [dt.date(2005, 1, 1), 4383],
"TAM06": [dt.date(2005, 1, 1), 4383],
"NVX": [dt.date(2005, 1, 1), 4383]}

# start date and number of daily time steps for gem lai file
site_dates_lai = {
"CAX04": [dt.date(2009, 1, 1), 730, 24],
"CAX06": [dt.date(2009, 1, 1), 730, 24],
"KEN01": [dt.date(2009, 1, 1), 1095, 36],
"KEN02": [dt.date(2009, 1, 1), 1095, 36],
"TAM05": [dt.date(2009, 1, 1), 1095, 36],
"TAM06": [dt.date(2009, 1, 1), 1461, 48]}

# start date and number of daily time steps for gem soil moisture data
site_years_sthuf = {
"KEN01": [2008, 2009, 2010, 2011, 2012, 2013],
}

if obs_data=="lba":
    path2csv="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/lba_obs/"
    path2nc="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/lba_obs/nc/"
elif obs_data=="gem_lai" or obs_data=="gem_sthuf":
    path2csv="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/GEM_Sophie/"
    path2nc="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/gem_obs/nc/"
elif obs_data=="gem_npp":
    path2nc="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/gem_obs/"	
    path2neemodel="/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/model_outputs/INLAND/site_runs/met/local/nc/"



# inland output variables
out_vars = ["swnet", "lwnet", "qle", "qh", "qg"] # EDIT THIS, add variable names, name them as you like

lba_obs_out_vars = ["NEE","NEEf","NEE_model","Re_5day_ust_Sco2_LUT","gpp_gb","par_fill","VPD","mrs"]
lba_obs_in_vars = ["NEE","NEEf","NEE_model","Re_5day_ust_Sco2_LUT","GEP_model","par_fill","VPD","mrs"]

# LAI, sthuf
gem_obs_out_vars = ["sthuf"]

if obs_data=="gem_npp":
	gem_obs_out_vars = ["root_NPP","total_litterfall_NPP","wood_NPP","branch_NPP","root_resp","resp_s_gb","wood_resp","Leaf_ flush","leaf_resp","root_NPP_se","total_litterfall_NPP_se","wood_NPP_se","branch_NPP_se","root_resp_se","resp_s_gb_se","wood_resp_se","Leaf_ flush_se","leaf_resp_se", "lai"]

#year,month,ingrowth_core_fine_root_NPP,total_litterfall_NPP,woody_NPP,branch_NPP,Rhizoshpere_respiration,Soil_heterotrophic_respiration,Wood_respiration,Leaf_ flush,Canopy_respiration,ingrowth_core_fine_root_NPP_se,total_litterfall_NPP_se,woody_NPP_se,branch_NPP_se,Rhizoshpere_respiration_se,Soil_heterotrophic_respiration_se,Wood_respiration_se,Leaf_ flush_se,Canopy_respiration_se

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

    carbon_lba_obs_df = pd.read_csv(filename, delimiter=',', header=None, names=CARBON_HEADERS)
    return carbon_lba_obs_df

    
def extract_npp_gem_obs(filename,site):

	print "Extracting NPP data for "+site
		
	NPP_HEADERS = ["year","month","ingrowth_core_fine_root_NPP","total_litterfall_NPP","woody_NPP","branch_NPP","Rhizoshpere_respiration","Soil_heterotrophic_respiration","Wood_respiration","Leaf_ flush","Canopy_respiration","ingrowth_core_fine_root_NPP_se","total_litterfall_NPP_se","woody_NPP_se","branch_NPP_se","Rhizoshpere_respiration_se","Soil_heterotrophic_respiration_se","Wood_respiration_se","Leaf_ flush_se","Canopy_respiration_se", "lai"]

	obs = pd.read_csv(filename)

	# re-set names of variables to shorter versions
	obs.rename(columns={"ingrowth_core_fine_root_NPP":"root_NPP"}, inplace=True)
	obs.rename(columns={"woody_NPP":"wood_NPP"}, inplace=True)
	obs.rename(columns={"Soil_heterotrophic_respiration":"resp_s_gb"}, inplace=True)
	obs.rename(columns={"Rhizoshpere_respiration":"root_resp"}, inplace=True)
	obs.rename(columns={"Wood_respiration":"wood_resp"}, inplace=True)
	obs.rename(columns={"Canopy_respiration":"leaf_resp"}, inplace=True)			

	# and their standard errors
	obs.rename(columns={"ingrowth_core_fine_root_NPP_se":"root_NPP_se"}, inplace=True)
	obs.rename(columns={"woody_NPP_se":"wood_NPP_se"}, inplace=True)
	obs.rename(columns={"Soil_heterotrophic_respiration_se":"resp_s_gb_se"}, inplace=True)
	obs.rename(columns={"Rhizoshpere_respiration_se":"root_resp_se"}, inplace=True)
	obs.rename(columns={"Wood_respiration_se":"wood_resp_se"}, inplace=True)
	obs.rename(columns={"Canopy_respiration_se":"leaf_resp_se"}, inplace=True)

	return obs



def extract_nee_obs(site):
	
	print "Extracting NEE for "+site

	# extract nee from obs
	fluxes = extract_carbon_lba_obs_csv("lba_obs/CAXday-carbon.dat",site)
	nee = fluxes["NEE_model"].tolist()
	nee_daily = pd.DataFrame({"NEE": nee}, columns=["NEE"])			

	nee_daily["date"] = pd.date_range(dt.date(1999,1,1), periods=731, freq="D") # add date
	nee_daily = nee_daily.set_index(["date"])
	
	nee_daily_updated = nee_daily.replace(-9999., np.nan) # replace -9999. with NaNs		
	nee_tmp = nee_daily_updated.resample('M').mean() # daily to monthly
	nee_tmp_list = nee_tmp["NEE"].tolist()
	
	# compute climatology
	nee_clim = []
	ni=0
	while ni<len(nee_tmp_list)/2:
		print ni
		nee_clim.append(np.nanmean([nee_tmp_list[ni],nee_tmp_list[ni+12]]))	
		ni+=1

	return nee_clim


def extract_nee_model(site):

	# extract nee from model
	ncfile = Dataset(path2neemodel+"CAX06-single_point-output_daily.nc",'r')
	nee = ncfile.variables['NEE'][:,0,0]

	model_tmp1 = pd.DataFrame(nee, columns=['NEE'])
	model_tmp1['date'] = pd.date_range(dt.date(2005,1,1), periods=len(nee), freq='D')
	model_tmp1 = model_tmp1.set_index(['date'])
	
	nee_tmp = model_tmp1.resample('M').mean() # daily to monthly	
	nee_tmp_list = nee_tmp["NEE"].tolist()	

	# compute climatology
	nee_clim = []
	ni=0
	while ni<len(nee_tmp_list)/12:
		print ni
		nee_clim.append(np.nanmean([nee_tmp_list[ni],nee_tmp_list[ni+12],nee_tmp_list[ni+24],nee_tmp_list[ni+36],nee_tmp_list[ni+48],nee_tmp_list[ni+60],nee_tmp_list[ni+72],nee_tmp_list[ni+84],nee_tmp_list[ni+96], nee_tmp_list[ni+108], nee_tmp_list[ni+120], nee_tmp_list[ni+132]]))	
		ni+=1

	return nee_clim





def write_nee_2netcdf(nee_obs,nee_model,site):

	# write to netCDF
	print "Writing LBA/GEM OBS output to netcdf with relevant metadata/attributes"

	data = Dataset(path2nc+site+"_nee.nc", 'w', format='NETCDF4')

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

	start_nee=dt.date(2005,1,1)
	times_all=[]
	for month in range(12):
		times_all.append(datetime.datetime.combine(start_nee, datetime.time())) 
		start_nee+=relativedelta(months=+1)	

	times=times_all
	time.units = 'days since 1850-01-01 00:00:00.0'
	time.calendar = 'gregorian'

	nee_obs_out=data.createVariable('nee_obs', 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)
	nee_model_out=data.createVariable('nee_inland', 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)

	nee_obs_out.units = "g.m-2.d-1"
	nee_obs_out.long_name = "Gridbox NEE"
	nee_obs_out.title = "NEE_GB from GEM"

	nee_obs_out[:,0,0] = nee_obs

	# min/max values
	nee_obs_out.actual_min = np.min(nee_obs_out[:,0,0])
	nee_obs_out.actual_max = np.max(nee_obs_out[:,0,0])

	nee_model_out.units = "g.m-2.d-1"
	nee_model_out.long_name = "Gridbox NEE"
	nee_model_out.title = "NEE_GB from INLAND"

	nee_model_out[:,0,0] = nee_model

	# min/max values
	nee_model_out.actual_min = np.min(nee_model_out[:,0,0])
	nee_model_out.actual_max = np.max(nee_model_out[:,0,0])

    	data.close()
    	print "***SUCCESS writing to file***"	
	







# tres is the temporal resolution of the output data
# LBA is hourly
# GEM can be daily or monthly
def write_obs_2netcdf(obs_output, site, outfilename, tres):
    """
    Takes a pandas dataframe with the obs data and writes it to netcdf, for a given site
    """

    # I modified this to work for both LBA/GEM data (DS)
    print "Writing LBA/GEM OBS output to netcdf with relevant metadata/attributes"

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

    time.calendar = 'gregorian'

    # add values to time variable
    if tres=="hourly":
    	times = [datetime.datetime.combine(site_dates[site][0], datetime.time()) + datetime.timedelta(hours=hour) for hour in range(site_dates[site][1]*24)]
	time.units = 'hours since 1850-01-01 00:00:00.0'

    elif (obs_data=="gem_lai" and tres=="monthly") or (obs_data=="gem_npp" and tres=="monthly"):

	if obs_data=="gem_npp":
		start_date=dt.date(2009,1,1)
		nyrs=2
	elif obs_data=="gem_lai":
		start_date=site_dates[site][0]
		nyrs=12

	times_all=[]
	for month in range(12*nyrs): # 12 years, 2005-2016, as I am dealing with GEM data
		times_all.append(datetime.datetime.combine(start_date, datetime.time())) 
		start_date+=relativedelta(months=+1)	

	times=times_all
	time.units = 'days since 1850-01-01 00:00:00.0'
	time.calendar = 'gregorian'
	
	
	#if obs_data=="gem_lai":
	#	data_tmp = obs_output["LAI"]
	#elif obs_data=="gem_npp":	
	#	data_tmp = obs_output["ingrowth_core_fine_root_NPP"]
		
	#data_updated = np.array([np.nan]*12*nyrs)

	#times_gem=[]
	
	#if obs_data=="gem_lai":
	#	start_date=site_dates_lai[site][0]
	#	for month in range(site_dates_lai[site][2]):
	#		times_gem.append(datetime.datetime.combine(start_date, datetime.time())) 
	#		start_date+=relativedelta(months=+1)
	#elif obs_data=="gem_npp":
	#	start_date=dt.date(2009,1,1)
	#	for month in range(731):
	#		times_gem.append(datetime.datetime.combine(start_date, datetime.time())) 
	#		start_date+=relativedelta(days=+1)

	

	# return indices of month data time object in daily list
	#month_inds=[]
	#i=0
	#while i<len(data):
	#	month_inds.append(times.index(times_gem[i]))		
	#	i+=1
		
	#if obs_data=="gem_lai":
	#	data_updated[month_inds]=data		
	#elif obs_data=="gem_npp":
	#	obs_output_updated={}
	#	for ob in obs_output:
	#		data_updated[month_inds]=data	

    
    elif obs_data=="gem_lai" and tres=="daily":

	times = [datetime.datetime.combine(site_dates[site][0], datetime.time()) + datetime.timedelta(days=day) for day in range(site_dates[site][1])]
	time.units = 'days since 1850-01-01 00:00:00.0'
	time.calendar = 'gregorian'

	# convert monthly LAI timeseries to daily
	# assume each LAI value was taken on the first of each month for now	
	lai = obs_output["LAI"]
	lai_updated = np.array([np.nan]*site_dates[site][1])

	# monthly time series of dates
	times_tmp=site_dates_lai[site][0]
	times_gem=[]
	for month in range(site_dates_lai[site][2]):
		times_gem.append(datetime.datetime.combine(times_tmp, datetime.time())) 
		times_tmp+=relativedelta(months=+1)

	# return indices of month data time object in daily list
	month_inds=[]
	i=0
	while i<len(lai):
		month_inds.append(times.index(times_gem[i]))		
		i+=1

	lai_updated[month_inds]=lai

    # saving time values
    print "Saving time values"	
    time[:] = date2num(times, time.units, calendar=time.calendar)

    print "Writing data to file for "+str(len(lats))+" point(s)"

    if obs_data=="lba":
    	obs_out_vars = lba_obs_in_vars
    elif obs_data=="gem_lai" or obs_data=="gem_npp":
	obs_out_vars = gem_obs_out_vars
    
    # calculate plant respiration ( root + wood + canopy respiration = autotrophic respiration)
    plant_resp = np.array(obs_output["root_resp"].tolist()) + np.array(obs_output["wood_resp"].tolist()) + np.array(obs_output["leaf_resp"].tolist())
    plant_resp_se = np.array(obs_output["root_resp_se"].tolist()) + np.array(obs_output["wood_resp_se"].tolist()) + np.array(obs_output["leaf_resp_se"].tolist())   

    dataout=data.createVariable("resp_p_gb", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)
    dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
    dataout.long_name = "plant respiration "	
    dataout.title = "plant respiration"

    dataout[:,0,0] = plant_resp
    
    dataout=data.createVariable("resp_p_gb_se", 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)
    dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
    dataout.long_name = "plant respiration se"	
    dataout.title = "plant respiration se"

    dataout[:,0,0] = plant_resp_se

    # min/max values
    dataout.actual_min = np.min(dataout[:,0,0])
    dataout.actual_max = np.max(dataout[:,0,0])

    for in_var in obs_out_vars:
	print "Creating netcdf variable: " + in_var
	
	if in_var=="LAI":
		dataout=data.createVariable('lai_gb', 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)
	else:
		dataout=data.createVariable(in_var, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)	


	# EDIT THIS, add metadata for other variables here
	if in_var=="GEP_model":  # GPP_GB alias	
		dataout.units = "gC.m-2.day-1"
		dataout.long_name = "Gridbox GPP"
		dataout.title = "GPP_GB"
	elif in_var=="lai":	
		dataout.units = "m2.m-2"
		dataout.long_name = "Gridbox leaf area index (LAI)"
		dataout.title = "LAI_GB"
	elif in_var=="root_NPP":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "ingrowth core fine root NPP"	
		dataout.title = "ingrowth core fine root NPP"
	elif in_var=="total_litterfall_NPP":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "total litterfall NPP. This leaf litterfall C and wood litterfall C."	
		dataout.title = "total litterfall NPP"
	elif in_var=="wood_NPP":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "woody NPP"	
		dataout.title = "woody NPP"
	elif in_var=="branch_NPP":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "branch NPP"	
		dataout.title = "branch NPP"
	elif in_var=="root_resp":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "Rhizoshpere_respiration"	
		dataout.title = "Rhizoshpere_respiration"
	elif in_var=="resp_s_gb":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "Soil heterotrophic respiration"	
		dataout.title = "Soil heterotrophic respiration"
	elif in_var=="wood_resp":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "Wood_respiration"	
		dataout.title = "Wood_respiration"
	elif in_var=="Leaf_flush":
		dataout.units = "MgC.ha-1.mo-1"
		dataout.long_name = "Leaf_flush"	
		dataout.title = "Leaf_flush"
	elif in_var=="leaf_resp":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "Canopy_respiration"	
		dataout.title = "Canopy_respiration"
	elif in_var=="root_NPP_se":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "ingrowth_core_fine_root_NPP_se"	
		dataout.title = "ingrowth_core_fine_root_NPP_se"			
	elif in_var=="total_litterfall_NPP_se":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "total_litterfall_NPP_se"	
		dataout.title = "total_litterfall_NPP_se"
	elif in_var=="wood_NPP_se":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "woody_NPP_se"	
		dataout.title = "woody_NPP_se"
	elif in_var=="branch_NPP_se":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "branch_NPP_se"	
		dataout.title = "branch_NPP_se"
	elif in_var=="root_resp_se":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "Rhizoshpere_respiration_se"	
		dataout.title = "Rhizoshpere_respiration_se"
	elif in_var=="resp_s_gb_se":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "Soil_heterotrophic_respiration_se"	
		dataout.title = "Soil_heterotrophic_respiration_se"
	elif in_var=="wood_resp_se":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "Wood_respiration_se"	
		dataout.title = "Wood_respiration_se"
	elif in_var=="Leaf_flush_se":
		dataout.units = "MgC.ha-1.mo-1"
		dataout.long_name = "Leaf_flush_se"	
		dataout.title = "Leaf_flush_se"
	elif in_var=="leaf_resp_se":
		dataout.units = "kg.m-2" #"MgC.ha-1.mo-1"
		dataout.long_name = "Canopy_respiration_se"	
		dataout.title = "Canopy_respiration_se"
	else:
		print "No data variables to write...? \n Check your netcdf file"

	print "Writing "+ in_var +" data to file..."
	
	if in_var=="root_NPP" or in_var=="total_litterfall_NPP" or in_var=="wood_NPP" or in_var=="branch_NPP" or in_var=="root_resp" or in_var=="resp_s_gb" or in_var=="wood_resp" or in_var=="leaf_resp" or in_var=="roott_NPP_se" or in_var=="total_litterfall_NPP_se" or in_var=="wood_NPP_se" or in_var=="branch_NPP_se" or in_var=="root_resp_se" or in_var=="resp_s_gb_se" or in_var=="wood_resp_se" or in_var=="leaf_resp_se":
		dataout[:,0,0] = obs_output[in_var].values*.1 # Mg.ha-1.month-1 -> kg.m-2			
	else:
		dataout[:,0,0] = obs_output[in_var].values

	# min/max values
	dataout.actual_min = np.min(dataout[:,0,0])
	dataout.actual_max = np.max(dataout[:,0,0])

    data.close()
    print "***SUCCESS writing to file***"	



def extract_gem_obs_csv(site):

	print "##############################################"
	print "Extracting GEM OBS data from csv file at " + site

	import pandas as pd

	# store data in a dictionary
	gem_data = {}

	for i in gem_obs_out_vars:

		if i=="LAI":
			gemfile = path2csv+"LAI/SouthAmerica_LAIfromCDoughty_July14_updated.csv"	
			print "Reading in "+gemfile

			# read in all data
			data = np.loadtxt(gemfile, delimiter=',', skiprows=1, dtype=np.str)

			plot_code = np.array(data[:,1])
			month = np.array(data[:,2])
			year = np.array(data[:,3])
			lai_all = np.array([float(x) for x in data[:,4]])

			# extract specific site data
			site_data = []
			j=0
			while j<len(plot_code):
				if site==plot_code[j]:
					site_data.append(lai_all[j])
				j+=1
			gem_data[i]=site_data

		elif i=="sthuf":
			gemfile = path2csv+"Soil_Moisture/all_plots.csv"

			print "Reading in "+gemfile
			
			# read in all data
			data = np.loadtxt(gemfile, delimiter=',', skiprows=1, dtype=np.str)

			plot_code = np.array(data[:,1])
			year = np.array([float(x) for x in data[:,2]])
			month = np.array([float(x) for x in data[:,3]])
			vwc = np.array([float(x) for x in data[:,5]])
			
			# extract data in 2 ways
			# first calculate a monthly average
			#site_data = {}	
			#j=0
			#while j<len(site_years_sthuf):
			#	k=0
			#	while k<len(range(1,12+1,1)):
			#			m=0
			#			while m<len(plot_code):							
			#				if site==plot_code[m] and :					

			#				m+=1
			#			k+=1
			#		k+=1		
			#	j+=1
			







	# create a pandas dataframe from a dictionary		
	gem_data_pd = pd.DataFrame.from_dict(gem_data)
	return gem_data_pd



##########################################################################################################
# extract data from all sites and write to netcdf
for site in sites:

	print "##############################################"
	print "PROCESSING SITE: " + site

	if obs_data=="lba":
    
		# Extract the data from the site csv
	    	#inland_outputs = extract_csv_data(path2csv+"/"+i+"-single_point-output.csv",i)
	    	lba_obs_outputs = extract_carbon_lba_obs_csv(path2csv+"/" + site + "day-carbon.dat", site)

	    	print "LBA OBS dataframe summary:"
	    	print lba_obs_outputs
	
	    	# Concatenate an output file csv
	    	outfile = path2nc+"/"+site+"-carbon_LBA_site_point.nc"

	    	# Write the variable to netcdf format
	    	write_obs_2netcdf(lba_obs_outputs, site, outfile,timeres)
	
	elif obs_data=="gem_lai":

		gem_obs_output = extract_gem_obs_csv(site)

	    	outfile = path2nc+"/"+site+"_GEM_site_point_"+timeres+".nc"

		print "Writing GEM data to netcdf"
		write_obs_2netcdf(gem_obs_output, site, outfile,timeres)

	elif obs_data=="gem_npp":
	
		gem_obs_output = extract_npp_gem_obs("Doughty_et_al_2015_CAX_Tower.csv",site)
	
		neeobs = extract_nee_obs(site)	
		neemodel = extract_nee_model(site)		

		write_nee_2netcdf(neeobs,neemodel,site)
		
		outfile = path2nc+"/"+site+"_GEM_"+timeres+".nc"

		print "Writing GEM NPP data to netcdf"
		write_obs_2netcdf(gem_obs_output, site, outfile,timeres)


##########################################################################################################


















