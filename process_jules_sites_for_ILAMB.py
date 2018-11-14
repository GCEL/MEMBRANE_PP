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
# GEM sites
# There is also a site called Tanguro
# Caxiuana: CAX04 (Control), CAX06 (Tower)
# Kenia: KEN01 (wet,deep), KEN02 (dry,shallow)
# Tambopata: TAM05, TAM06
# Nova Xavantina: NVX
GEM_sites=["CAX04", "CAX06", "KEN01", "KEN02", "TAM05", "TAM06", "NVX"]

# LBA sites
# BAN: Bananal Island site
# FNS: Fazenda Nossa Senhora
# K34: Manaus k34
# K67: Santarem k67
# K77: Santarem k77
# K83: Santarem k83
# PDG: Reserva Pe-de-Gigante
# RJA: Reserva Jaru Forest
LBA_sites=["BAN", "FNS", "K34", "K67", "K77", "K83", "PDG", "RJA"]

##########################################################################################################
# add sites to this array for model runs
sites=["BAN"]
#sites=["BAN", "FNS", "K34", "K67", "K77", "K83", "PDG", "RJA", "CAX04", "CAX06", "KEN01", "KEN02", "TAM05", "TAM06", "NVX"]

# TRMM, MSWEP
# only for the GEM sites as there is no local precip data available
precip_data="MSWEP"

# Met data can be local (LBA, GEM) or global (CRUJRA)
met_data="global"

# model simulation type
# MET-LOCAL: LBA or GEM met data is used to drive the model.
# MET-GLOBAL: CRUJRA met data is used to drive the model.
# LAI: Model runs in which lai is prescribed. Local data used.
# STHUF: Model runs in which soil moisture (sthuf) is prescribed. Local data used. 
sim_type="MET-GLOBAL"

# specifies output folders
output_type={"MET-GLOBAL": "CRUJRA", "LAI": "PRESC_LAI", "STHUF": "PRESC_STHUF"}

# model output temporal resolution
# global: 6hourly, daily, monthly, yearly
# local: hourly, daily, yearly
output_timesteps=["daily"]

output_tres={"hourly": "H", "3hourly": "3H", "6hourly": "6H", "daily": "D", "monthly": "M", "yearly": "Y"}

# number of soil layers
nlayers = 4

path2inputs="/exports/csce/datastore/geos/users/dslevin/jules/5.2/"

##########################################################################################################
# site: [lat, lon]
sitelocs = {"CAX04": [-1.716, -51.457], "CAX06": [-1.737, -51.462], "KEN01": [-16.016, -62.730], "KEN02": [-16.016, -62.730],
"TAM05": [-12.831, -69.271], "TAM06": [-12.839, -69.296], "NVX": [-14.70, -52.35], "BAN": [-9.82, -50.15], "FNS": [-10.76, -62.36],
"K34": [-2.61, -60.21], "K67": [-2.85, -54.97], "K77": [-3.01, -54.54], "K83": [-3.05, -54.93],
"PDG": [-21.62, -47.63], "RJA": [-10.08, -62.36]}

# start date and number of daily time steps
site_dates = {"BAN": [dt.date(2004, 1, 2), 1033], "FNS": [dt.date(1999, 1, 2), 1094], "K34": [dt.date(2003, 1, 2), 1017], "K67": [dt.date(2002, 1, 2), 685],\
"K77": [dt.date(2001, 1, 2), 1824], "K83": [dt.date(2001, 1, 2), 952], "PDG": [dt.date(2002, 1, 2), 728], "RJA": [dt.date(2000, 2, 3), 953],\
"CAX04": [dt.date(2005, 1, 1), 4383], "CAX06": [dt.date(2005, 1, 1), 4383], "KEN01": [dt.date(2005, 1, 1), 4383], "KEN02": [dt.date(2005, 1, 1), 4383],\
"TAM05": [dt.date(2005, 1, 1), 4383], "TAM06": [dt.date(2005, 1, 1), 4383], "NVX": [dt.date(2005, 1, 1), 4383]}

# frac (Fractional cover of each surface type), 
# gpp_gb (gridbox GPP), resp_p_gb (gridbox plant respiration), resp_s_gb (gridbox soil respiration), 
# npp_n_gb (NPP (GBM) post N-limitation), cs_gb (gridbox soil carbon, total), 
# cv (Gridbox mean vegetation carbon at end of model timestep)
# lai (PFT leaf area index),
# leafC (PFT carbon in leaf biomass (kg m-2 )), woodC (PFT carbon in woody biomass (kg m-2 )), 
# rootC (PFT carbon in root biomass (kg m-2 )),
# lit_c_mean (Gridbox mean carbon litter (kg m-2 (360days)-1))
# leaf_litC (Litter CARBON due to leaf turnover (kg/m2/360 days))
# root_litC (Litter CARBON due to root turnover (kg/m2/360 days))
# wood_litC (Litter CARBON due to wood turnover (kg/m2/360 days))

# EDIT THIS, jules output variables, you have to use these names as you are reading the variables from netcdf files
if sim_type=="MET-LOCAL" or sim_type=="MET-GLOBAL":
	varlist = ["gpp_gb", "resp_p_gb", "resp_s_gb", "npp_n_gb", "lai", "leafC", "woodC", "rootC", "lit_c_mean", "leaf_litC", "root_litC", "wood_litC", "latent_heat", "fsmc_gb", "smc_avail_top", "fqw_gb", "et_stom_gb", "sthu", "sthf"]
elif sim_type=="LAI":
	varlist = ["gpp_gb", "resp_p_gb", "resp_s_gb", "lai", "latent_heat", "fsmc_gb", "smc_avail_top", "fqw_gb", "et_stom_gb", "sthu", "sthf"]
elif sim_type=="STHUF":
	varlist = ["gpp_gb", "resp_p_gb", "resp_s_gb", "lai", "latent_heat", "fsmc_gb", "smc_avail_top", "fqw_gb", "et_stom_gb", "sthu", "sthf"]

##########################################################################################################
# convert PFT values to gridbox average
def convert_pft_2gridbox(pft_data,fracs):

	grid_data=[]

	i=0
	while i<len(pft_data[:,0]):

		grid_data.append(np.sum(pft_data[i,:]*fracs[i,:]))

		i+=1

	return np.array(grid_data)



def extract_data(filename):

	print "Extracting JULES output from "+filename
	
	data = Dataset(filename, "r")

	jules_data = {}

	i=0
	while i<len(varlist):

		###########################################################################################
		# EDIT THIS	
		###########################################################################################
		if varlist[i]=="gpp_gb" or varlist[i]=="resp_p_gb" or varlist[i]=="resp_s_gb" or varlist[i]=="npp_gb":
			jules_data[varlist[i]] = np.array(data.variables[varlist[i]][:][:,0,0])*1000*86400 # kgC.m-2.s-1 -> gC.m-2.day-1
		elif varlist[i]=="cs_gb" or varlist[i]=="cv" or varlist[i]=="npp_n_gb":
			jules_data[varlist[i]] = np.array(data.variables[varlist[i]][:][:,0,0]) # kgC.m-2	
		elif varlist[i]=="lai":
			pft_fracs = np.array(data.variables["frac"][:][:,:,0,0]) # dimensionless (0-1)	
			lai_pft = np.array(data.variables[varlist[i]][:][:,:,0,0]) # m2.m-2

			# convert PFT lai to gridbox lai
			jules_data[varlist[i]] = convert_pft_2gridbox(lai_pft,pft_fracs[:,0:9])
		elif varlist[i]=="leafC" or varlist[i]=="woodC" or varlist[i]=="rootC":
			pft_fracs = np.array(data.variables["frac"][:][:,:,0,0]) # dimensionless (0-1)	
			biomass_pft = np.array(data.variables[varlist[i]][:][:,:,0,0]) # kgC.m-2

			# convert PFT biomass to gridbox biomass
			jules_data[varlist[i]] = convert_pft_2gridbox(biomass_pft,pft_fracs[:,0:9])
		elif varlist[i]=="lit_c_mean":
			jules_data[varlist[i]] = np.array(data.variables[varlist[i]][:][:,0,0]) # kgC.m-2.360days-1
		elif varlist[i]=="leaf_litC" or varlist[i]=="root_litC" or varlist[i]=="wood_litC":
			pft_fracs = np.array(data.variables["frac"][:][:,:,0,0]) # dimensionless (0-1)	
			biomass_pft = np.array(data.variables[varlist[i]][:][:,:,0,0]) # kgC.m-2

			# convert to gridbox values
			jules_data[varlist[i]] = convert_pft_2gridbox(biomass_pft,pft_fracs[:,0:9])
		elif varlist[i]=="latent_heat" or varlist[i]=="fsmc_gb" or varlist[i]=="smc_avail_top" or varlist[i]=="fqw_gb" or varlist[i]=="et_stom_gb":
			jules_data[varlist[i]] = np.array(data.variables[varlist[i]][:][:,0,0])
		elif varlist[i]=="sthu" or varlist[i]=="sthf":
			jules_data[varlist[i]] = np.array(data.variables[varlist[i]][:][:,:,0,0])



		i+=1		

	data.close()
	
	return jules_data 


def write2netCDF(moutput, site, outfilename, tres):

	print "Writing JULES output (with metadata, evaluation tool ready) to file"
	
	data = Dataset(outfilename, 'w', format='NETCDF4')
	
	lats = [sitelocs[site][0]]
	lons = [sitelocs[site][1]]
	
	# dimensions
	latdim=data.createDimension('lat', len(lats))
	londim=data.createDimension('lon', len(lons))
	tdim=data.createDimension('time', None) # record, or unlimited dimension
	sdim=data.createDimension('soil', nlayers)

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

	# extract overlapping data from runs using local and global met data
	if tres=="daily" and sim_type=="MET-GLOBAL":



	elif tres=="hourly" and :
		# add values to time variable
		times = [datetime.datetime.combine(site_dates[site][0], datetime.time()) + datetime.timedelta(hours=hour) for hour in range(site_dates[site][1]*24)]
 	 
	    	time.units = 'hours since 1850-01-01 00:00:00.0'
	    	time.calendar = 'gregorian'

	time[:] = date2num(times, time.units, calendar=time.calendar)

	print "Writing data to file for "+str(len(lats))+" point(s)"

	for i in varlist:

		# EDIT THIS, add metadata for other variables here
		if i=="sthu" or i=="sthf": # these variables have outputs for each soil layer
			dataout=data.createVariable(i, 'f4', ('time', 'soil', 'lat', 'lon',),fill_value=np.nan)
		else:
			dataout=data.createVariable(i, 'f4', ('time', 'lat', 'lon',),fill_value=np.nan)

		if i=="gpp_gb":	
			dataout.units = "gC.m-2.day-1"
			dataout.long_name = "Gridbox GPP"
		elif i=="resp_p_gb":	
			dataout.units = "gC.m-2.day-1"
			dataout.long_name = "Gridbox Plant respiration"
		elif i=="resp_s_gb":	
			dataout.units = "gC.m-2.day-1"
			dataout.long_name = "Gridbox Total soil respiration"
		elif i=="npp_n_gb":
			dataout.units = "gC.m-2.day-1"
			dataout.long_name = "NPP (GBM) post N-limitation"
		elif i=="cs_gb":
			dataout.units = "kgC.m-2"
			dataout.long_name = "gridbox soil carbon, total"
		elif i=="cv":
			dataout.units = "kgC.m-2"
			dataout.long_name = "Gridbox mean vegetation carbon at end of model timestep"			
		elif i=="lai":
			dataout.units = "m2.m-2"
			dataout.long_name = "Gridbox leaf area index"
		elif i=="leafC":
			dataout.units = "KgC.m-2"
			dataout.long_name = "Gridbox carbon in leaf biomass"
		elif i=="woodC":
			dataout.units = "KgC.m-2"
			dataout.long_name = "Gridbox carbon in woody biomass"
		elif i=="rootC":		
			dataout.units = "KgC.m-2"
			dataout.long_name = "Gridbox carbon in root biomass"
		elif i=="lit_c_mean":
			dataout.units = "Kg.m-2.360days"
			dataout.long_name = "Gridbox mean carbon litter"
		elif i=="leaf_litC":
			dataout.units = "KgC.m-2.360days"
			dataout.long_name = "Litter CARBON due to leaf turnover"
		elif i=="root_litC":
			dataout.units = "KgC.m-2.360days"
			dataout.long_name = "Litter CARBON due to root turnover"
		elif i=="wood_litC":
			dataout.units = "KgC.m-2.360days"
			dataout.long_name = "Litter CARBON due to wood turnover"
		elif i=="latent_heat":
			dataout.units = "W.m-2"
			dataout.long_name = "Gridbox surface latent heat flux"
		elif i=="fsmc_gb":
			dataout.units = "dimensionless (0-1)"
			dataout.long_name = "Gridbox soil moisture availability factor (beta)"
		elif i=="smc_avail_top":
			dataout.units = "Kg.m-2"
			dataout.long_name = "Gridbox available moisture in top 1.00000000m of soil"
		elif i=="fqw_gb":
			dataout.units = "Kg.m-2.s-1"
			dataout.long_name = "Gridbox moisture flux from surface"
		elif i=="et_stom_gb":
			dataout.units = "Kg.m-2.s-1"
			dataout.long_name = "Gridbox stomatal transpiration"
		elif i=="sthu":
			dataout.units = "dimensionless (0-1)"
			dataout.long_name = "Gridbox unfrozen moisture content of each soil layer as a fraction of saturation"
		elif i=="sthf":
			dataout.units = "dimensionless (0-1)"
			dataout.long_name = "Gridbox frozen moisture content of each soil layer as a fraction of saturation"

		print "Writing "+i+" data to file..."

		if i=="sthu" or i=="sthf":
			dataout[:,:,0,0] = moutput[i][:,:]
			
			# min/max values
			dataout.actual_min = np.min(dataout[:,:,0,0])
			dataout.actual_max = np.max(dataout[:,:,0,0])			
		else:
			dataout[:,0,0] = moutput[i][:]
		
			# min/max values
			dataout.actual_min = np.min(dataout[:,0,0])
			dataout.actual_max = np.max(dataout[:,0,0])

	data.close()
	print "***SUCCESS writing to file***"	


##########################################################################################################
print "Post-processing "+sim_type+" simulations"
i=0
while i<len(sites):

	print "####################################################################################"

	# select model simulation type
	if sim_type=="MET-LOCAL":
		if sites[i] in LBA_sites:
			
			for tstep in output_timesteps:				

				print "Processing "+sites[i]+" ["+sim_type+", "+tstep+"]"

				infile = path2inputs+"jules_vn5.2_doc/examples/"+sites[i]+"_9PFTS/output/"+sites[i]+"_LBA_DYN."+output_tres[tstep]+".nc"
				outfile = path2inputs+"jules_vn5.2_doc/examples/"+sites[i]+"_9PFTS/output/"+sites[i]+"_LBA_DYN."+output_tres[tstep]+".pp.nc"

				model_output = extract_data(infile)
				write2netCDF(model_output,sites[i],outfile,tstep)

		elif sites[i] in GEM_sites:
	
			for tstep in output_timesteps:

				print "Processing "+sites[i]+" ["+sim_type+", "+tstep+"]"

				infile = path2inputs+"jules_vn5.2_doc/examples/"+sites[i]+"_"+precip_data+"_9PFTS/output/"+sites[i]+"_GEM_DYN."+output_tres[tstep]+".nc"
				outfile = path2inputs+"jules_vn5.2_doc/examples/"+sites[i]+"_"+precip_data+"_9PFTS/output/"+sites[i]+"_GEM_DYN."+output_tres[tstep]+".pp.nc"

				model_output = extract_data(infile)
				write2netCDF(model_output,sites[i],outfile,tstep)

	elif sim_type=="MET-GLOBAL":
			
		for tstep in output_timesteps:

			print "Processing "+sites[i]+" ["+sim_type+", "+tstep+"]"

			infile = path2inputs+"jules_vn5.2_doc/examples/"+sites[i]+"_"+output_type[sim_type]+"_9PFTS/output/"+sites[i]+"_CRUJRA_DYN."+output_tres[tstep]+".nc"
			outfile = path2inputs+"jules_vn5.2_doc/examples/"+sites[i]+"_"+output_type[sim_type]+"_9PFTS/output/"+sites[i]+"_CRUJRA_DYN."+output_tres[tstep]+".pp.nc"		

			model_output = extract_data(infile)
			write2netCDF(model_output,sites[i],outfile,tstep)
	
	elif sim_type=="LAI" or sim_type=="STHUF":
		
		if sites[i] in LBA_sites:

			for tstep in output_timesteps:

				print "Processing "+sites[i]+" ["+sim_type+", "+tstep+"]"

				infile = path2inputs+"jules_vn5.2_doc/examples/"+sites[i]+"_"+output_type[sim_type]+"_9PFTS/output/"+sites[i]+"_LBA_DYN."+output_tres[tstep]+".nc"
				outfile = path2inputs+"jules_vn5.2_doc/examples/"+sites[i]+"_"+output_type[sim_type]+"_9PFTS/output/"+sites[i]+"_LBA_DYN."+output_tres[tstep]+".pp.nc"

				model_output = extract_data(infile)
				write2netCDF(model_output,sites[i],outfile,tstep)

		elif sites[i] in GEM_sites:

			for tstep in output_timesteps:

				print "Processing "+sites[i]+" ["+sim_type+", "+tstep+"]"

				infile = path2inputs+"jules_vn5.2_doc/examples/"+sites[i]+"_"+output_type[sim_type]+"_9PFTS/output/"+sites[i]+"_GEM_DYN."+output_tres[tstep]+".nc"
				outfile = path2inputs+"jules_vn5.2_doc/examples/"+sites[i]+"_"+output_type[sim_type]+"_9PFTS/output/"+sites[i]+"_GEM_DYN."+output_tres[tstep]+".pp.nc"

				model_output = extract_data(infile)
				write2netCDF(model_output,sites[i],outfile,tstep)		
	

		







	i+=1





##########################################################################################################


















     
