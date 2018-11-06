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
# LBA sites
# BAN: Bananal Island site
# FNS: Fazenda Nossa Senhora
# K34: Manaus k34
# K67: Santarem k67
# K77: Santarem k77
# K83: Santarem k83
# PDG: Reserva Pe-de-Gigante
# RJA: Reserva Jaru Forest

# GPP, Reco, NEE, H, LE
julesvar = "NEE"

julesvar_model = {"GPP": "gpp_gb", "Reco_plant": "resp_p_gb", "Reco_soil": "resp_s_gb"}

sites = ["BAN", "FNS", "K34", "K67", "K77", "K83", "PDG", "RJA"]

site_years = {"BAN": range(2004,2006+1,1), "FNS": range(1999,2001+1,1), "K34": range(2003,2005+1,1), "K67": range(2002,2003+1,1),\
"K77": range(2001,2005+1,1), "K83": range(2001,2003+1,1), "PDG": range(2002,2003+1,1), "RJA": range(2000,2002+1,1)}

# start date and number of daily time steps
site_dates = {"BAN": [dt.date(2004, 1, 2), 1033], "FNS": [dt.date(1999, 1, 2), 1094], "K34": [dt.date(2003, 1, 2), 1017], "K67": [dt.date(2002, 1, 2), 685], "K77": [dt.date(2001, 1, 2), 1824], "K83": [dt.date(2001, 1, 2), 952], "PDG": [dt.date(2002, 1, 2), 728], "RJA": [dt.date(2000, 2, 3), 953]}

##########################################################################################################
def extract_model_output():

	print "Extracting JULES output from file"

	jules_data = {}
	
	i=0
	while i<len(sites):

		print "Extracting data for "+sites[i]

		filename = "jules_vn5.2_doc/examples/"+sites[i]+"_9PFTS/output/LBA_"+sites[i]+".D.nc"
		data = Dataset(filename, "r")

		if julesvar=="GPP":
			jules_data[sites[i]] = np.array(data.variables["gpp_gb"][:][:,0,0])*1000*86400 # kgC.m-2.s-1 -> gC.m-2.day-1
		elif julesvar=="Reco":	
			resp_p_gb = np.array(data.variables["resp_p_gb"][:][:,0,0])*1000*86400
			resp_s_gb = np.array(data.variables["resp_s_gb"][:][:,0,0])*1000*86400		
			jules_data[sites[i]] = resp_p_gb + resp_s_gb # reco 
		elif julesvar=="NEE":
			gpp = np.array(data.variables["gpp_gb"][:][:,0,0])*1000*86400
			resp_p = np.array(data.variables["resp_p_gb"][:][:,0,0])*1000*86400
			resp_s = np.array(data.variables["resp_s_gb"][:][:,0,0])*1000*86400
			jules_data[sites[i]] = (resp_p + resp_s) - gpp # nee			
	
		data.close()
		i+=1

	return jules_data 	




# create daily timeseries of dates using start date and number of days
def create_dates_daily(start,num_days):

	# keep a copy of start date
	start_first = start

	dates_daily = [] # store daily dates

	j=0
	while j<num_days:

		delta_m = relativedelta(days=1)
		start += delta_m
		dates_daily.append(start)

		j+=1

	# remove last and add first
	#dates_daily.pop()
	
	return dates_daily


# read in obs
def read_obs(var):

	print "Reading in LBA obs"

	obs_lba = {}; times = {}

	i=0
	while i<len(sites):

		filename = "jules_vn5.2_doc/data_for_suite_u-al752/vn0.5/lba_obs/"+sites[i]+"day-carbon.dat"
		data = np.loadtxt(filename, delimiter=',', skiprows=1, dtype='str')	
	
		if var=="GPP":
			obs_lba_tmp = data[:,7]
		elif var=="Reco":
			obs_lba_tmp = data[:,6]
		elif var=="NEE":
			obs_lba_tmp = data[:,5]

		obs_lba_tmp = np.array([float(x) for x in obs_lba_tmp])
		obs_lba_tmp[obs_lba_tmp==-9999.]=np.nan

		if var in ["GPP", "Reco", "NEE"]:
			obs_lba[sites[i]] = obs_lba_tmp*1.0E-6*12*86400 # g.m-2.day

		# read in times as well
		years = data[:,0]
		
		times[sites[i]] = create_dates_daily(dt.date(int(years[0]), 1, 1),len(years))
	
		i+=1
	return obs_lba, times  








##########################################################################################################
jules_output = extract_model_output()

lba_obs, lba_times = read_obs(julesvar)

##########################################################################################################
print "Plotting data"
fig=plt.figure(0,figsize=(20,10)) # figsize = width x height
adj = plt.subplots_adjust(hspace=0.5,wspace=0.2)

# line type: solid line, dashed line, symbols, etc.
ltype_solid = "-"
ltype_dashed = "--"
ltype_symbols = "^"
ltype_dotted = ":"

model_c = '#B22400'
obs_c = '#077FD8'

# Plot data for each variable
i=0
while i<len(sites):

	print "Plotting "+julesvar+" data"

	print "Plotting "+sites[i]
	ax=fig.add_subplot(len(sites)/2,2,i+1)

	years = YearLocator() # every year
	months = MonthLocator() # every month
	yearsFmt = DateFormatter('%Y')

	sim_dates = create_dates_daily(site_dates[sites[i]][0], site_dates[sites[i]][1])

	# use only obs which overlaps with model output
	if lba_times[sites[i]][0] < sim_dates[0]:
		start_ind = 0
	else:
		start_ind = sim_dates.index(lba_times[sites[i]][0])
			

	if lba_times[sites[i]][len(lba_times[sites[i]])-1] > sim_dates[len(sim_dates)-1]: 	
		end_ind = len(sim_dates)-1
	else:
		end_ind = sim_dates.index(lba_times[sites[i]][len(lba_times[sites[i]])-1])
	
	jules_output_updated = jules_output[sites[i]][start_ind:end_ind+1]
	lba_obs_updated = lba_obs[sites[i]][start_ind:end_ind+1]

	datemin = min(sim_dates)
	datemax = max(sim_dates)

	# set min/max for x-axis
	ax.set_xlim(datemin, datemax)

	print len(sim_dates), len(jules_output_updated)
	print len(sim_dates), len(lba_obs_updated)

	ax.plot_date(sim_dates,jules_output_updated,ltype_solid,linewidth=1.5,color=model_c,label="JULES5.2")	
	ax.plot_date(sim_dates,lba_obs_updated,ltype_solid,linewidth=1.5,color=obs_c,label="obs")

	ax.set_title(sites[i], fontsize=18)

	# add axes labels, limits and title
	if julesvar=="GPP":
		if i==0 or i==2 or i==4 or i==6:
			plt.ylabel('$\mathregular{gC\,m^{-2}day^{-1}}$',fontsize=14)
		ax.set_ylim(0,20)	
	elif julesvar=="Reco":
		if i==0 or i==2 or i==4 or i==6:
			plt.ylabel('$\mathregular{gC\,m^{-2}day^{-1}}$',fontsize=14)
		ax.set_ylim(0,20)
	elif julesvar=="NEE":
		if i==0 or i==2 or i==4 or i==6:
			plt.ylabel('$\mathregular{gC\,m^{-2}day^{-1}}$',fontsize=14)
		ax.set_ylim(-10,10)
		ax.axhline(linewidth=1.2, color='k')

	
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
fig.savefig("plots/"+julesvar+"_LBA.png", format="png", dpi=300)
plt.show()











































