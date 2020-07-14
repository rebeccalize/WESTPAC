import sys
import fnmatch
import os
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,cm
import matplotlib.cm as cm
import numpy as np
from matplotlib.patches import Polygon
from netCDF4 import Dataset, num2date, netcdftime
import datetime
import matplotlib.path as mpath
import matplotlib
from scipy.interpolate import spline

#This script produces the mean daily rainfall, and total daily rainfall, using TRMM, for each phase of the MJO (amplitude >= 1)
#The ONI stuff at the start is, I think, when I started to look at ENSO-MJO combinations, but changed the methodology and also started to use IMERG instead of TRMM
#So this stuff isn't actually used... 


TRMM_DIR = "/gws/nopw/j04/klingaman/datasets/TRMM_3B42/V7_NC_daily/"

ONIfile = "/home/users/emerton/nino3.4.ONI.NCAR-climate-data-guide.3-month_running_mean.1979-2020.txt"
ONI_data = np.genfromtxt(ONIfile, dtype=float, skip_header=8)

print ONI_data[:-1,0]
print len(ONI_data[:-1,0])

ONI_ts = []

#move the data from a row per year with columns going from DJF to NDJ, to a continuous timeseries of 3-month-averages from 1979 to 2020
for i in range(len(ONI_data[:-1,0])):
	#print ONI_data[i,1:]
	for j in range(1,13):
		ONI_ts.append(ONI_data[i,j])
	
#print ONI_ts

#np.concatenate(ONI_ts)

ONI_EN = list(ONI_ts)
ONI_LN = list(ONI_ts)
ONI_n = list(ONI_ts)

#split the data up by value so that we can colour the line to show el nino and la nina events
for i in range(len(ONI_ts)-1):
	if ONI_ts[i] >= 0.5:
		ONI_EN[i-1] = ONI_ts[i-1]
		ONI_EN[i] = ONI_ts[i]
		ONI_EN[i+1] = ONI_ts[i+1]
	else:
		ONI_EN[i] = np.nan
		
for i in range(len(ONI_ts)-1):
	if ONI_ts[i] <= -0.5:
		ONI_LN[i] = ONI_ts[i]
		ONI_LN[i-1] = ONI_ts[i-1]
		ONI_LN[i+1] = ONI_ts[i+1]
		
	else:
		ONI_LN[i] = np.nan
		
for i in range(len(ONI_ts)-1):	
	if -0.5 <= ONI_ts[i] <= 0.5:
		ONI_n[i] = ONI_ts[i]
		ONI_n[i-1] = ONI_ts[i-1]
		ONI_n[i+1] = ONI_ts[i-1]
		
	else:
		ONI_n[i] = np.nan
		
#print ONI_EN
#print ONI_LN
#print ONI_n

years = [1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018] #

	
for MJO in [8]: #1,2,3,4,5,6,7,8
	MJO_file = "/gws/nopw/j04/klingaman/datasets/MJO_INDICES/MJO_phase"+str(MJO)+".jan-dec_dmeans_ts.1979-2019.nc"
	ffMJO = Dataset(MJO_file,'r')
	MJOampdata = ffMJO.variables['phase_ts'][:]
	MJOdatesnc = ffMJO.variables['time'][:]

	t_unit = ffMJO.variables['time'].units
	t_cal = ffMJO.variables['time'].calendar
	tvalue = num2date(MJOdatesnc,units=t_unit, calendar=t_cal)
	
	MJOdates = [i.strftime("%Y-%m-%d") for i in tvalue]
	
	#print MJOdates
	

	trmm_comp = np.zeros((400,1440)) #trmm grid - use this array to compute the average precip per day in this MJO phase
	trmm_sum = np.zeros((400,1440)) #trmm grid - use this array to hold the sum of the daily precip on days in this MJO phase
	no_days = 0 #number of days included in thesum, in order to compute the mean 

	for year in years: #
	
		trmm_file = TRMM_DIR+"3B42_daily."+str(year)+".nc"
		print trmm_file
		
		ff_trmm = Dataset(trmm_file,'r')
		
		if year >= 2016: #variable name changes from 'r' to 'precipitation' from 2016 onwards, and has no time variable so need to create it
			daily_rain = ff_trmm.variables['precipitation'][:]
			
			start_time = datetime.datetime(year,1,1)
			trmm_tvalue = np.array([start_time + datetime.timedelta(days=i) for i in xrange(len(daily_rain[:,0,0]))])
			
			trmm_dates = [j.strftime("%Y-%m-%d") for j in trmm_tvalue]
			
		else: #for pre-2016, precip variable is called r, and we can use the existing time variable to get the dates
			daily_rain = ff_trmm.variables['r'][:]
		
			trmmdatesnc = ff_trmm.variables['time'][:]
			trmm_t_unit = ff_trmm.variables['time'].units
			trmm_t_cal = u"standard"
		
			trmm_tvalue = num2date(trmmdatesnc, units = trmm_t_unit, calendar=trmm_t_cal)
		
			trmm_dates = [j.strftime("%Y-%m-%d") for j in trmm_tvalue]
		
		print trmm_dates
		
		for date in trmm_dates:
			
			dateinfo = datetime.datetime.strptime(date, "%Y-%m-%d")
			if dateinfo.strftime("%m-%d") == '02-29': #leap year 29th Febs missing in MJO data
			
				continue
				
			else:
				i = MJOdates.index(date)  #np.where(MJOdates = date)
		
				mjo_amp = MJOampdata[i]
			
				print mjo_amp
			
				if mjo_amp >= 1.0:
			
					a = trmm_dates.index(date)
				
					print a
			
					trmm_sum[:,:] += daily_rain[a,:,:]
				
					trmm_sum = np.add(trmm_sum, daily_rain[a,:,:])
				
					no_days += 1
				
				else:
			
					continue
				
			
	print trmm_sum
	print no_days
	
	for x in range(400):
		for y in range(1440):
			trmm_comp[x,y] = trmm_sum[x,y] / no_days
			
			
	np.savetxt("trmm_precip_composite_mean_daily_rainfall.MJO_phase_"+str(MJO)+"."+str(years[0])+"-"+str(years[-1])+"."+str(no_days).zfill(4)+"-days.txt",  trmm_comp[:,:], '%.4f')
	np.savetxt("trmm_precip_composite_sum_of_daily_rainfall.MJO_phase_"+str(MJO)+"."+str(years[0])+"-"+str(years[-1])+"."+str(no_days).zfill(4)+"-days.txt", trmm_sum[:,:], '%.4f')
			
	
			
			
		



