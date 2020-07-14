import sys
#sys.path.append('/usr/lib/python2.7/site-packages/mpl_toolkits/')
#sys.path.append('/usr/lib/python2.7/site-packages/')
import fnmatch
import os
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap, cm
#import mpl_toolkits
#mpl_toolkits.__path__.append('/usr/lib/python2.7/site-packages/mpl_toolkits/')
from mpl_toolkits.basemap import Basemap,cm
import matplotlib.cm as cm
import numpy as np
from matplotlib.patches import Polygon
from netCDF4 import Dataset, num2date, netcdftime
import datetime
import matplotlib.path as mpath
import matplotlib
from scipy.interpolate import spline


TRMM_DIR = "/gws/nopw/j04/klingaman/datasets/TRMM_3B42/V7_NC_daily/"


TRMM_CLIM_DIR = "/gws/nopw/j04/klingaman/emerton/"
TRMM_CLIM_FILE = TRMM_CLIM_DIR+"TRMM_3B42.daily_climatology.2000-2018.nc"

print TRMM_CLIM_FILE

ffCLIM = Dataset(TRMM_CLIM_FILE,'r')

DAILY_CLIM = ffCLIM.variables['r'][:]
			
start_time_clim = datetime.datetime(2000,1,1)
trmm_tvalue_clim = np.array([start_time_clim + datetime.timedelta(days=i) for i in xrange(len(DAILY_CLIM[:,0,0]))])
			
trmm_dates_clim = [j.strftime("%Y-%m-%d") for j in trmm_tvalue_clim]




years = [2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018] #,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018

	
for MJO in [4]: #1,2,3,4,5,6,7,8
	MJO_file = "/gws/nopw/j04/klingaman/datasets/MJO_INDICES/MJO_phase"+str(MJO)+".jan-dec_dmeans_ts.1979-2019.nc"
	ffMJO = Dataset(MJO_file,'r')
	MJOampdata = ffMJO.variables['phase_ts'][:]
	MJOdatesnc = ffMJO.variables['time'][:]

	t_unit = ffMJO.variables['time'].units
	t_cal = ffMJO.variables['time'].calendar
	tvalue = num2date(MJOdatesnc,units=t_unit, calendar=t_cal)
	
	MJOdates = [i.strftime("%Y-%m-%d") for i in tvalue]
	
	#print MJOdates
	

	trmm_bias_comp = np.zeros((400,1440)) #trmm grid - use this array to compute the average precip per day in this MJO phase
	trmm_bias_array = np.zeros((400,1440,1)) #use vstack with this to hold the bias values, as don't want to sum; trmm grid - use this array to hold the sum of the daily precip on days in this MJO phase
	
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
			
				#print mjo_amp
			
				if mjo_amp >= 1.0:
				
					no_days += 1
				
					print date
			
					a = trmm_dates.index(date)
					
					#the climatology file uses default dates from the start year of 2000, so get the month and day of this date but use 2000 to find the index of this date in the climatology data
					b = trmm_dates_clim.index("2000-"+dateinfo.strftime("%m")+"-"+dateinfo.strftime("%d"))
					
					bias = np.subtract(daily_rain[a,:,:], DAILY_CLIM[b,:,:])
					
					print "bias shape: ", np.shape(bias)
					
					if no_days == 1:
					
						trmm_bias_array[:,:,0] = bias[:,:]
							
					else:
						#dstack stacks up 2d arrays using the third dimension, i.e. stacking up images
						trmm_bias_array = np.dstack((trmm_bias_array, bias))
						
	
					
				
				else:
			
					continue
				

	print np.shape(trmm_bias_array)
	print no_days
	
	for x in range(400):
		for y in range(1440):
		
			trmm_bias_comp[x,y] = np.nanmean(trmm_bias_array[x,y,:])
			
			
	np.savetxt("trmm_precip_composite.daily_rainfall_bias_from_daily_climatology.MJO_phase_"+str(MJO)+"."+str(years[0])+"-"+str(years[-1])+"."+str(no_days).zfill(4)+"-days.txt",  trmm_bias_comp[:,:], '%.4f')
	#np.savetxt("trmm_precip_composite_sum_of_daily_rainfall.MJO_phase_"+str(MJO)+"."+str(years[0])+"-"+str(years[-1])+"."+str(no_days).zfill(4)+"-days.txt", trmm_sum[:,:], '%.4f')
			
	
			
			
		



