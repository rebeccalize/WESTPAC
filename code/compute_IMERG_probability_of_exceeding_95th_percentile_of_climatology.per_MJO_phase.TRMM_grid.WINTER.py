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

#This version uses the 95th percentile that is calculated across the whole year (rather than the 95th percentile for a given day across all years)
#And all values below 1mm of rainfall were removed before calculating the 95th percentile

IMERG_DIR = "/gws/nopw/j04/klingaman/emerton/"


IMERG_DIR = "/gws/nopw/j04/klingaman/emerton/GPM-IMERG/"

MONTHLY_CLIM_ARRAY = np.zeros((6,1440,400))

IMERG_CLIM_DIR = "/gws/nopw/j04/klingaman/emerton/GPM-IMERG/"

for m,i in zip([10,11,12,1,2,3],[0,1,2,3,4,5]):


	IMERG_CLIM_FILE = IMERG_CLIM_DIR+"GPM_IMERG.DAILY.TRMM_GRID.2001-2019.m"+str(m).zfill(2)+".95th_percentile.nc"


	print IMERG_CLIM_FILE

	ffCLIM = Dataset(IMERG_CLIM_FILE,'r')

	CLIM = ffCLIM.variables['precipitationCal'][:]
	
	MONTHLY_CLIM_ARRAY[i,:,:] = CLIM

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
		print np.shape(MONTHLY_CLIM_ARRAY)
			

		winter_start_date = imerg_dates.index(str(year)+'-10-01')
	
	
		#we want to run this only for dates in winter - so we want to run October 1st to March 31st
		#This first section does October 1st to December 31st
		#Then below, we do Jan 1st to March 31st, as it crosses the end of the file so we can't just loop through them
		
		#Found the index for the first of August, loop over this date to the end of the array	
		for date in imerg_dates[winter_start_date:len(daily_rain[:,0,0])]:
		
		
			#print date
			
			dateinfo = datetime.datetime.strptime(date, "%Y-%m-%d")
			
			if dateinfo.strftime("%m") == '10':
				ci = 0
			elif dateinfo.strftime("%m") == '11':
				ci = 1
			elif dateinfo.strftime("%m") == '12':
				ci = 2
			elif dateinfo.strftime("%m") == '01':
				ci = 3
			elif dateinfo.strftime("%m") == '02':
				ci = 4
			elif dateinfo.strftime("%m") == '03':
				ci = 5
			else:
				print "month is ", dateinfo.strftime("%m")
			
			if dateinfo.strftime("%m-%d") == '02-29': #leap year 29th Febs missing in MJO data
			
				continue
				
			else:
				i = MJOdates.index(date)  #np.where(MJOdates = date)
		
				mjo_amp = MJOampdata[i]
			
				#print "mjo amp: ",mjo_amp
			
				if mjo_amp >= 1.0:
				
					total_no_days += 1
				
					print date
			
					a = imerg_dates.index(date)
					
					for x in range(1440):
						for y in range(400):
							
							#if the daily rainfall on this day and at this location,is equal to or above the 95th percentile of climatology for this day, 
							#add one to the array that holds the total number of days exceedingthe 95th percentile at each gridpoint
							
							if daily_rain[a,x,y] >= MONTHLY_CLIM_ARRAY[ci,x,y]:
								no_days_exceeding_pctl[x,y] += 1
								
							#otherwise, do nothing, we're only interested in finding out how many days exceeded the 95th percentile during this MJO phase
							else:
								continue
							
								
			
						
	
					
				#if the MJO amplitude for this phase is <1 on this day, we're not interested so carry on
				else:
			
					continue
				

		winter_end_date = imerg_dates.index(str(year)+'-03-31')
	
		#Found the index for the last day of July, loop over the start of the dates until then, i.e. first half of the year	
		for date in imerg_dates[0:winter_end_date]:
		
			#print date
			
			dateinfo = datetime.datetime.strptime(date, "%Y-%m-%d")
			
			if dateinfo.strftime("%m") == '10':
				ci = 0
			elif dateinfo.strftime("%m") == '11':
				ci = 1
			elif dateinfo.strftime("%m") == '12':
				ci = 2
			elif dateinfo.strftime("%m") == '01':
				ci = 3
			elif dateinfo.strftime("%m") == '02':
				ci = 4
			elif dateinfo.strftime("%m") == '03':
				ci = 5
			else:
				print "month is ", dateinfo.strftime("%m")
			
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
					
					for x in range(1440):
						for y in range(400):
							
							#if the daily rainfall on this day and at this location,is equal to or above the 95th percentile of climatology for this day, 
							#add one to the array that holds the total number of days exceedingthe 95th percentile at each gridpoint
							
							if daily_rain[a,x,y] >= MONTHLY_CLIM_ARRAY[ci,x,y]:
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
			
	np.savetxt("gpm-imerg.precip_composite.percentage_of_days_exceeding_95th_percentile_of_WINTER_climatology.MJO_phase_"+str(MJO)+"."+str(years[0])+"-"+str(years[-1])+"."+str(total_no_days).zfill(4)+"-days.using_precipitationCal.lowerres_TRMM_grid.Oct-March.MONTHLY_percentiles.txt",  perc_days_exc_pctl[:,:], '%.4f')
	
						
			
			
	
	
			
	
			
			
		



