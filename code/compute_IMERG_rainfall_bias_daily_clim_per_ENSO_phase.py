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


IMERG_DIR = "/gws/nopw/j04/klingaman/emerton/"


IMERG_CLIM_DIR = "/gws/nopw/j04/klingaman/emerton/"
IMERG_CLIM_FILE = IMERG_CLIM_DIR+"GPM_IMERG.daily_climatology.2001-2019_regridded_to_TRMM_grid.nc"

print IMERG_CLIM_FILE

ffCLIM = Dataset(IMERG_CLIM_FILE,'r')

DAILY_CLIM = ffCLIM.variables['precipitationCal'][:] #precipitationCal
			
start_time_clim = datetime.datetime(2001,1,1)
imerg_tvalue_clim = np.array([start_time_clim + datetime.timedelta(days=i) for i in xrange(len(DAILY_CLIM[:,0,0]))])
			
imerg_dates_clim = [j.strftime("%Y-%m-%d") for j in imerg_tvalue_clim]

ffCLIM.close()


#Older version, using threshold of 0.5 for the ONI (Nick suggested to limit to >= 1.0 to avoid very weak events such as 2014/2015)
#ENyears = [2002, 2004, 2006, 2009, 2014, 2015, 2018] #start year of the El Nino, take months from August of this year through to July of the next
#LNyears = [2005, 2007, 2008, 2010, 2011, 2016, 2017]
#Nyears = [2001, 2003, 2012, 2013]

#Newer version, using a threshold of 1.0 for the ONI, to only include somewhat stronger events (weaker threshold includes e.g. 2014/2015 "bust")
ENyears=[2002, 2009, 2015]
LNyears=[2007, 2010, 2011, 2017]
Nyears=[2001,2003,2004,2005,2006,2008,2012,2013,2014,2016,2018]

	
def compute_rainfall_anomalies(year_array): #1,2,3,4,5,6,7,8

	
	imerg_bias_comp = np.zeros((1440,400)) #imerg grid - use this array to compute the average precip per day in this MJO phase
	imerg_bias_array = np.zeros((1440,400,1)) #use vstack with this to hold the bias values, as don't want to sum; imerg grid - use this array to hold the sum of the daily precip on days in this MJO phase
	
	no_days = 0 #number of days included in thesum, in order to compute the mean 

	for year in year_array: #
	
		imerg_file = IMERG_DIR+"GPM_IMERG_3B-DAY.MS.MRG.3IMERG."+str(year)+".DAILY.V06_regridded_to_TRMM_grid.nc"
		
		print imerg_file
		
		ff_imerg = Dataset(imerg_file,'r')
		
		daily_rain = ff_imerg.variables['precipitationCal'][:] #precipitationCal
			
		start_time = datetime.datetime(year,1,1)
		imerg_tvalue = np.array([start_time + datetime.timedelta(days=i) for i in xrange(len(daily_rain[:,0,0]))])
			
		imerg_dates = [j.strftime("%Y-%m-%d") for j in imerg_tvalue]
		
		ff_imerg.close()
		
		print imerg_dates
		
		#El Nino / La Nina events tend to peak in winter, developing from summer, peaking in winter and decaying into the next spring/summer
		#So no point looking at a Jan - December year
		#Instead, we want to look at August - December of this year, and January to July of the following year
		a = imerg_dates.index(str(year)+'-08-01')
	
		#Found the index for the first of August, loop over this date to the end of the array	
		for date in imerg_dates[a:len(daily_rain[:,0,0])]:
		
		
			print date
			
			print no_days
			
			dateinfo = datetime.datetime.strptime(date, "%Y-%m-%d")
			
			if dateinfo.strftime("%m-%d") == '02-29': #leap year 29th Febs missing in MJO data
			
				continue
				
			else:
				no_days += 1 
			
				a = imerg_dates.index(date)
					
				#the climatology file uses default dates from the start year of 2001, so get the month and day of this date but use 2000 to find the index of this date in the climatology data
				b = imerg_dates_clim.index("2001-"+dateinfo.strftime("%m")+"-"+dateinfo.strftime("%d"))
					
				bias = np.subtract(daily_rain[a,:,:], DAILY_CLIM[b,:,:])
					
				print "bias shape: ", np.shape(bias)
					
				if no_days == 1:
				
					imerg_bias_array[:,:,0] = bias[:,:]
							
				else:
					#dstack stacks up 2d arrays using the third dimension, i.e. stacking up images
					imerg_bias_array = np.dstack((imerg_bias_array, bias))
						
					
					
		#aaaaaaaand repeat for the first half of the following year!
					
					
		imerg_file_2 = IMERG_DIR+"GPM_IMERG_3B-DAY.MS.MRG.3IMERG."+str(year+1)+".DAILY.V06_regridded_to_TRMM_grid.nc"
		
		print imerg_file_2
		
		ff_imerg_2 = Dataset(imerg_file_2,'r')
		
		daily_rain_2 = ff_imerg_2.variables['precipitationCal'][:] #precipitationCal
			
		start_time_2 = datetime.datetime(year+1,1,1)
		imerg_tvalue_2 = np.array([start_time_2 + datetime.timedelta(days=i) for i in xrange(len(daily_rain_2[:,0,0]))])
			
		imerg_dates_2 = [j.strftime("%Y-%m-%d") for j in imerg_tvalue_2]
		
		ff_imerg_2.close()
		
		print imerg_dates_2
		
		#El Nino / La Nina events tend to peak in winter, developing from summer, peaking in winter and decaying into the next spring/summer
		#So no point looking at a Jan - December year
		#Instead, we want to look at August - December of this year, and January to July of the following year
		b = imerg_dates_2.index(str(year+1)+'-07-31')
	
		#Found the index for the last day of July, loop over the start of the dates until then, i.e. first half of the year	
		for date in imerg_dates_2[0:b]:
		
			print date
			
			print no_days
			
			dateinfo = datetime.datetime.strptime(date, "%Y-%m-%d")
			
			if dateinfo.strftime("%m-%d") == '02-29': #leap year 29th Febs missing in MJO data
			
				continue
				
			else:
			
				no_days += 1
			
				a = imerg_dates_2.index(date)
					
				#the climatology file uses default dates from the start year of 2001, so get the month and day of this date but use 2000 to find the index of this date in the climatology data
				b = imerg_dates_clim.index("2001-"+dateinfo.strftime("%m")+"-"+dateinfo.strftime("%d"))
					
				bias = np.subtract(daily_rain_2[a,:,:], DAILY_CLIM[b,:,:])
					
				print "bias shape: ", np.shape(bias)
					
				if no_days == 1:
					
					imerg_bias_array[:,:,0] = bias[:,:]
							
				else:
					#dstack stacks up 2d arrays using the third dimension, i.e. stacking up images
					imerg_bias_array = np.dstack((imerg_bias_array, bias))
						
				

	print np.shape(imerg_bias_array)
	print no_days

	
	for x in range(1440):
		for y in range(400):
		
			imerg_bias_comp[x,y] = np.nanmean(imerg_bias_array[x,y,:])
			
	if year_array==ENyears:
		enso_label='elnino_years_1.0_ONI_threshold'
	elif year_array==LNyears:
		enso_label='lanina_years_1.0_ONI_threshold'
			
			
	np.savetxt("gpm-imerg.precip_composite.daily_rainfall_bias_from_daily_climatology."+enso_label+"."+str(no_days).zfill(4)+"-days.txt",  imerg_bias_comp[:,:], '%.4f')
	
	
	#in this last bit, count how many days out of all the MJO days have the same sign anomaly as the mean
	#then convert the number of days with the same sign as the mean, to a % of the total number of days, and save this
			
	consistency=np.zeros((1440,400))
	
	for i in range(len(imerg_bias_array[0,0,:])):
		
		for x in range(1440):
			for y in range(400):
				
				if imerg_bias_comp[x,y] == 0:
					
					if imerg_bias_array[x,y,i] == 0:
					
						consistency[x,y] += 1
						
					else:
						continue
						
				elif imerg_bias_comp[x,y] > 0:
					
					if imerg_bias_array[x,y,i] > 0:
					
						consistency[x,y] += 1
					else:
						continue
						
				elif imerg_bias_comp[x,y] < 0:
				
					if imerg_bias_array[x,y,i] < 0:
					
						consistency[x,y] += 1
						
					else:
						continue
						
						
	for x in range(1440):
		for y in range(400):
			
			consistency[x,y] = (consistency[x,y]/no_days)*100
			
	
	print imerg_bias_comp
	print np.nanmean(imerg_bias_comp)
	print np.nanmax(imerg_bias_comp)
	print np.nanmean(imerg_bias_comp)
	
	

		
		
	
	np.savetxt("gpm-imerg.CONSISTENCY.precip_composite.daily_rainfall_bias_from_daily_climatology."+enso_label+"."+str(no_days).zfill(4)+"-days.txt",  consistency[:,:], '%.4f')
	
			
			
compute_rainfall_anomalies(ENyears)
#compute_rainfall_anomalies(LNyears)


