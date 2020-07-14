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


#For days in each MJO phase, count the number of days on which the rainfall exceeded the 95th percentile of climatology for that day
#Convert this to the percentage of days on which the rainfall exceeded the 95th percentile of climatology for that day

#We can then later plot the bias in probability i.e. the probability in this MJO phase - 5% (the climatological probability of exceeding the 95th percentile is 5%)
#To see whether the probability is increased or decreased during each phase of the MJO, of having extreme rainfall 


IMERG_DIR = "/gws/nopw/j04/klingaman/emerton/"


IMERG_CLIM_DIR = "/gws/nopw/j04/klingaman/emerton/"
IMERG_CLIM_FILE = IMERG_CLIM_DIR+"GPM_IMERG.daily_extremevals_climatology.2001-2019.95th_percentile.TRMM_grid.nc"

print IMERG_CLIM_FILE

ffCLIM = Dataset(IMERG_CLIM_FILE,'r')

DAILY_CLIM = ffCLIM.variables['precipitationCal'][:]
			
start_time_clim = datetime.datetime(2001,1,1)
imerg_tvalue_clim = np.array([start_time_clim + datetime.timedelta(days=i) for i in xrange(len(DAILY_CLIM[:,0,0]))])
			
imerg_dates_clim = [j.strftime("%Y-%m-%d") for j in imerg_tvalue_clim]

ffCLIM.close()



years = [2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019] #,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018

	
for MJO in [8]: #1,2,3,4,5,6,7,8
	MJO_file = "/gws/nopw/j04/klingaman/datasets/MJO_INDICES/MJO_phase"+str(MJO)+".jan-dec_dmeans_ts.1979-2019.nc"
	ffMJO = Dataset(MJO_file,'r')
	MJOampdata = ffMJO.variables['phase_ts'][:]
	MJOdatesnc = ffMJO.variables['time'][:]

	t_unit = ffMJO.variables['time'].units
	t_cal = ffMJO.variables['time'].calendar
	tvalue = num2date(MJOdatesnc,units=t_unit, calendar=t_cal)
	
	MJOdates = [i.strftime("%Y-%m-%d") for i in tvalue]
	
	ffMJO.close()
	
	#print MJOdates
	
	#array to hold the total number of days in this MJO phase, at each gridpoint, that had rainfall exceeding the 95th percentile of climatology
	no_days_exceeding_pctl = np.zeros((1440,400))
	
	#we also want to know the total number of days in this MJO phase, so we can calculate what % of them had extreme rainfall
	total_no_days = 0 #number of days included in thesum, in order to compute the mean 

	for year in years: #
	
		imerg_file = IMERG_DIR+"GPM_IMERG_3B-DAY.MS.MRG.3IMERG."+str(year)+".DAILY.V06_regridded_to_TRMM_grid.nc"
		print imerg_file
		
		ff_imerg = Dataset(imerg_file,'r')
		
		daily_rain = ff_imerg.variables['precipitationCal'][:]
			
		start_time = datetime.datetime(year,1,1)
		imerg_tvalue = np.array([start_time + datetime.timedelta(days=i) for i in xrange(len(daily_rain[:,0,0]))])
			
		imerg_dates = [j.strftime("%Y-%m-%d") for j in imerg_tvalue]
		
		ff_imerg.close()
		
		print np.shape(daily_rain)
		print np.shape(DAILY_CLIM)
			

		
		print imerg_dates
		
		for date in imerg_dates:
			
			dateinfo = datetime.datetime.strptime(date, "%Y-%m-%d")
			
			if dateinfo.strftime("%m-%d") == '02-29': #leap year 29th Febs missing in MJO data
			
				continue
				
			else:
				i = MJOdates.index(date)  #np.where(MJOdates = date)
		
				mjo_amp = MJOampdata[i]
			
				#print mjo_amp
			
				if mjo_amp >= 1.0:
				
					total_no_days += 1
				
					print date
			
					a = imerg_dates.index(date)
					
					#the climatology file uses default dates from the start year of 2001, so get the month and day of this date but use 2000 to find the index of this date in the climatology data
					b = imerg_dates_clim.index("2001-"+dateinfo.strftime("%m")+"-"+dateinfo.strftime("%d"))
					
					
					#GREAT! we have found a day where the MJO amplitude is greater than 1 for this phase,
					#and we've found the index of that date in the rainfall data for this year, and in the 95th percentile of climatology file
					#now, does the rainfall on this day exceed the 95th percentile of climatology for this day?
					#at each point...
					
					#bias = np.subtract(daily_rain[a,:,:], DAILY_CLIM[b,:,:])
					
					for x in range(1440):
						for y in range(400):
							
							#if the daily rainfall on this day and at this location,is equal to or above the 95th percentile of climatology for this day, 
							#add one to the array that holds the total number of days exceedingthe 95th percentile at each gridpoint
							
							if daily_rain[a,x,y] >= DAILY_CLIM[b,x,y]:
								no_days_exceeding_pctl[x,y] += 1
								
							#otherwise, do nothing, we're only interested in finding out how many days exceeded the 95th percentile during this MJO phase
							else:
								continue
							
								
			
						
	
					
				#if the MJO amplitude for this phase is <1 on this day, we're not interested so carry on
				else:
			
					continue
				

	
	print total_no_days
	print no_days_exceeding_pctl
	print np.nanmax(no_days_exceeding_pctl)
	
	
	#convert the number of days exceeding the 95th percnetile, to the % of days exceeding the 95th percentils
	
	perc_days_exc_pctl = np.zeros((1440,400))
	
	for x in range(1440):
		for y in range(400):
		
			perc_days_exc_pctl[x,y] = (no_days_exceeding_pctl[x,y] / total_no_days) * 100
			
	np.savetxt("gpm-imerg.precip_composite.percentage_of_days_exceeding_95th_percentile_of_climatology.MJO_phase_"+str(MJO)+"."+str(years[0])+"-"+str(years[-1])+"."+str(total_no_days).zfill(4)+"-days.using_precipitationCal_lowerres_TRMM_grid.txt",  perc_days_exc_pctl[:,:], '%.4f')
	
						
			
			
	
	
			
	
			
			
		



