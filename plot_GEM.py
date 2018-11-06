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

from netCDF4 import Dataset

import calendar
import datetime as dt
from dateutil.relativedelta import relativedelta

##########################################################################################################
# Caxiuana, Kenia, Nova Xavantina, Tambopata
# Caxiuana: CAX04, CAX06
# Kenia: KEN01, KEN02
# Tambopata:  TAM05, TAM06
# Nova Xavantina:  NVX
julesvar = "gpp_gb"

sites = ["CAX04", "KEN01", "TAM05", "NVX"]

precipdata = "TRMM"

years = range(2005,2016+1,1)

# day (D)
tres="D"

# gpp_gb = gridbox gpp
# resp_p_gb = gridbox plant respiration
# resp_s_gb = gridbox total soil respiration 
# resp_l = PFT leaf respiration
# resp_w = PFT wood respiration
# resp_r = PFT root respiration

# c_fluxes, water
vartype = "c_fluxes"

varlist_C_fluxes = ["gpp_gb", "resp_p_gb", "resp_s_gb", "resp_l", "resp_w", "resp_r"]

varlist_water = ["smc_tot", "fsmc_gb"]

if vartype=="c_fluxes":
	varlist = varlist_C_fluxes
elif vartype=="water":
	varlist = varlist_water

##########################################################################################################
def extract_data(filename):

	print "Extracting JULES output from file"
	
	data = Dataset(filename, "r")

	if julesvar=="gpp_gb" or julesvar=="resp_p_gb" or julesvar=="resp_s_gb":
		jules_data = np.array(data.variables[julesvar][:][:,0,0])*1000*86400 # kgC.m-2.s-1 -> gC.m-2.day-1
	elif julesvar=="resp_l" or julesvar=="resp_w" or julesvar=="resp_r":		
		jules_data = np.array(data.variables[julesvar][:][:,0,0,0])*1000*86400 # kgC.m-2.s-1 -> gC.m-2.day-1
	elif julesvar=="latent_heat" or julesvar=="fsmc_gb" or julesvar=="smc_tot" or julesvar=="fsmc_gb" or julesvar=="smc_avail_top":
		jules_data = np.array(data.variables[julesvar][:][:,0,0])

	data.close()
	
	return jules_data 	


# create daily timeseries of dates
def create_dates_daily(start):

	# keep a copy of start date
	start_first = start

	dates_daily = [] # store daily dates

	num_days = 0.	
	for i in years:
		if calendar.isleap(i):
			num_days+=366
		else:
			num_days+=365 

	j=0
	while j<num_days:

		delta_m = relativedelta(days=1)
		start += delta_m
		dates_daily.append(start)

		j+=1

	# remove last and add first
	dates_daily.pop()
	
	return dates_daily

##########################################################################################################
# Extract data for all sites
site_data={}

i=0
while i<len(sites):

	print "Extracting data for "+sites[i]
	site_data[sites[i]]=extract_data("jules_vn5.2_doc/examples/"+sites[i]+"_"+precipdata+"_9pfts/output/"+sites[i]+"_"+precipdata+"_9pfts."+tres+".nc")

	i+=1

##########################################################################################################
print "Plotting data"
fig=plt.figure(0,figsize=(20,10)) # figsize = width x height
adj = plt.subplots_adjust(hspace=0.5,wspace=0.2)

d0 = dt.date(years[0], 1, 1)

model_dates = create_dates_daily(d0)

# line type: solid line, dashed line, symbols, etc.
ltype_solid = "-"
ltype_dashed = "--"
ltype_symbols = "^"
ltype_dotted = ":"

model_c = '#B22400'

# Plot data for each variable
i=0
while i<len(sites):

	print "Plotting "+julesvar+" data"

	print "Plotting "+sites[i]
	ax=fig.add_subplot(len(sites),1,i+1)

	years = YearLocator() # every year
	months = MonthLocator() # every month
	yearsFmt = DateFormatter('%Y')

	ax.plot(model_dates,site_data[sites[i]],ltype_solid,linewidth=1.5,color=model_c,label="JULES5.2")	

	ax.set_title(sites[i], fontsize=18)

	# add axes labels, limits and title
	if julesvar=="gpp_gb":
		plt.ylabel('$\mathregular{gC\,m^{-2}day^{-1}}$',fontsize=14)
		ax.set_ylim(0,20)	
	elif julesvar=="resp_p_gb":
		plt.ylabel('$\mathregular{gC\,m^{-2}day^{-1}}$',fontsize=14)
		ax.set_ylim(0,12)
	elif julesvar=="resp_s_gb":
		plt.ylabel('$\mathregular{gC\,m^{-2}day^{-1}}$',fontsize=14)
		ax.set_ylim(0,15)
	elif julesvar=="resp_l":
		plt.ylabel('$\mathregular{gC\,m^{-2}day^{-1}}$',fontsize=14)
		ax.set_ylim(0,7)
		ax.set_title("PFT leaf respiration", fontsize=18)
	elif julesvar=="resp_w":
		plt.ylabel('$\mathregular{gC\,m^{-2}day^{-1}}$',fontsize=14)
		ax.set_ylim(0,7)
		ax.set_title("PFT wood respiration", fontsize=18)
	elif julesvar=="resp_r":
		plt.ylabel('$\mathregular{gC\,m^{-2}day^{-1}}$',fontsize=14)
		ax.set_ylim(0,7)
		ax.set_title("PFT root respiration", fontsize=18)
	elif julesvar=="latent_heat":
		plt.ylabel('$\mathregular{Wm^{-2}}$',fontsize=14)
		ax.set_ylim(0,400)
		ax.set_title("Latent heat", fontsize=18)
	elif julesvar=="fsmc_gb":
		plt.ylabel('Beta factor',fontsize=14)
		#ax.set_ylim(0,1)
		ax.set_title("Gridbox soil moisture availability factor (beta)", fontsize=18)
	elif julesvar=="smc_tot":
		plt.ylabel('$\mathregular{kg\,m^{-2}}$',fontsize=14)
		#ax.set_ylim(0,10)
		ax.set_title("Gridbox total soil moisture in column", fontsize=18)


	# format the ticks
	ax.xaxis.set_major_locator(years)
	ax.xaxis.set_major_formatter(yearsFmt)
	ax.xaxis.set_minor_locator(months)

	# set min/max for x-axis
	#ax.set_xlim(datemin, datemax)

	for tick in ax.xaxis.get_minor_ticks():
	    tick.tick1line.set_markersize(4)
	    tick.tick2line.set_markersize(0)
	    tick.label1.set_horizontalalignment('center')
	    tick.label.set_fontsize(12)

	# Change fontsize of x-axis labels
	for tick in ax.xaxis.get_major_ticks():
	  tick.label.set_fontsize(11)

	# Change fontsize of y-axis labels
	for tick in ax.yaxis.get_major_ticks():
	  tick.label.set_fontsize(11)

	# add legends
	if i==0:

		# labelspacing = vertical spacing between legend labels
		# columnspacing = spacing between columns 
		legend = ax.legend(bbox_to_anchor=(1., 1.08), numpoints=1, handletextpad=0.5, frameon=False, labelspacing=.05, columnspacing=.5, ncol=3)

		# Set fontsize
		for label in legend.get_texts():
    			label.set_fontsize(14)

	i+=1

# Add global title
fig.suptitle(julesvar, fontsize=24)

##########################################################################################################
print "Saving figure as png"
fig.savefig("plots/"+julesvar+"_GEM.png", format="png", dpi=300)
plt.show()











































