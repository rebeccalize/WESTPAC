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

DAILY_CLIM = ffCLIM.variables['precipitationCal'][:]
			
start_time_clim = datetime.datetime(2001,1,1)
imerg_tvalue_clim = np.array([start_time_clim + datetime.timedelta(days=i) for i in xrange(len(DAILY_CLIM[:,0,0]))])
			
imerg_dates_clim = [j.strftime("%Y-%m-%d") for j in imerg_tvalue_clim]

ffCLIM.close()


BSISO_FILE = "/gws/nopw/j04/klingaman/datasets/BSISO_INDICES/bsiso_indices.jan-dec_dmeans_ts.1981-2019.nc"

print BSISO_FILE

ffBSISO = Dataset(BSISO_FILE,'r')
BSISOdatesnc = ffBSISO.variables['time'][:]

BSISO_amp_data = ffBSISO.variables['bsiso1_amp'][:]
BSISO_phase_data = ffBSISO.variables['bsiso1_phase'][:]


#I want the phase and amplitude of PAIR 1 (?!) of the BSISO, as there are two BSISO indices? see Lee et al., 2013
#The file has some variables I'm unsure of, i.e. 'Component 1 / Component 2 of the BSISO 1 index" - I think these come together to give the overall amplitude that I use, like the MJO
#and also amp_0 of BSISO 1 index. What is the difference between amp and amp_0? The variable descriptions are the same. I have used 'amp' for now. 

t_unit = ffBSISO.variables['time'].units
t_cal = ffBSISO.variables['time'].calendar
tvalue = num2date(BSISOdatesnc, units=t_unit, calendar=t_cal)

BSISOdates = [i.strftime("%Y-%m-%d") for i in tvalue]

ffBSISO.close()

years = [2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019] #,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018

#So for the BSISO, I don't have separate file for each phase. 
#So first we need to check "if BSISO_phase = phase_of_interest" then "if amplitude > 1", do stuff
	
for BSISO in [6]: #1,2,3,4,5,6,7,8

	imerg_bias_comp = np.zeros((1440,400)) #imerg grid - use this array to compute the average precip per day in this MJO phase
	imerg_bias_array = np.zeros((1440,400,1)) #use vstack with this to hold the bias values, as don't want to sum; imerg grid - use this array to hold the sum of the daily precip on days in this MJO phase
	
	no_days = 0 #number of days included in thesum, in order to compute the mean 

	for year in years: #
	
		imerg_file = IMERG_DIR+"GPM_IMERG_3B-DAY.MS.MRG.3IMERG."+str(year)+".DAILY.V06_regridded_to_TRMM_grid.nc"
		print imerg_file
		
		ff_imerg = Dataset(imerg_file,'r')
		
		daily_rain = ff_imerg.variables['precipitationCal'][:]
			
		start_time = datetime.datetime(year,1,1)
		imerg_tvalue = np.array([start_time + datetime.timedelta(days=i) for i in xrange(len(daily_rain[:,0,0]))])
			
		imerg_dates = [j.strftime("%Y-%m-%d") for j in imerg_tvalue]
		
		ff_imerg.close()
			

		
		print imerg_dates
		
		summer_start_date = imerg_dates.index(str(year)+'-04-01')
		summer_end_date = imerg_dates.index(str(year)+'-09-30')
	
	
		#we want to run this only for dates in 'extended summer' - so we want to run April 1st to September 30th
		#For winter, we needed two loops since this crosses the end of the file, but don't need to do this for summer
		
		#Found the index for the first of August, loop over this date to the end of the array	
		for date in imerg_dates[summer_start_date:summer_end_date]:
		
			#print date
			
			dateinfo = datetime.datetime.strptime(date, "%Y-%m-%d")
			
			if dateinfo.strftime("%m-%d") == '02-29': #leap year 29th Febs missing in BSISO data
			
				continue
				
			else:
				i = BSISOdates.index(date)  #np.where(MJOdates = date)
		
				BSISO_phase = BSISO_phase_data[i]
				BSISO_amp = BSISO_amp_data[i]
			
				#print mjo_amp
				
				if BSISO_phase == BSISO and BSISO_amp >= 1.0:
				
					print BSISO_phase
					print BSISO_amp
				
					no_days += 1
				
					print date
			
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
						
	
					
				
				else:
			
					continue
					
		

	print np.shape(imerg_bias_array)
	print no_days
	print len(imerg_bias_array[0,0,:])
	#the above two should be the same, or the consistency calcs won't work?
	
	
	consistency=np.zeros((1440,400))
	
	for x in range(1440):
		for y in range(400):
		
			imerg_bias_comp[x,y] = np.nanmean(imerg_bias_array[x,y,:])
			
	np.savetxt("gpm-imerg.SUMMER_April-Sept.precip_composite.daily_rainfall_bias_from_daily_climatology.BSISO_phase_"+str(BSISO)+"."+str(years[0])+"-"+str(years[-1])+"."+str(no_days).zfill(4)+"-days.lowerres_TRMM_grid.txt",  imerg_bias_comp[:,:], '%.4f')
	
	#in this last bit, count how many days out of all the MJO days have the same sign anomaly as the mean
	#then convert the number of days with the same sign as the mean, to a % of the total number of days, and save this
			
	for i in range(len(imerg_bias_array[0,0,:])):
		
		for x in range(1440):
			for y in range(400):
				
				if imerg_bias_comp[x,y] == 0:
					
					if imerg_bias_array[x,y,i] == 0:
					
						consistency[x,y] += 1
						
					else:
						continue
						
				if imerg_bias_comp[x,y] > 0:
					
					if imerg_bias_array[x,y,i] > 0:
					
						consistency[x,y] += 1
					else:
						continue
						
				if imerg_bias_comp[x,y] < 0:
				
					if imerg_bias_array[x,y,i] < 0:
					
						consistency[x,y] += 1
						
					else:
						continue
						
						
	for x in range(1440):
		for y in range(400):
			
			consistency[x,y] = (consistency[x,y]/no_days)*100
						
			
			
	
	np.savetxt("gpm-imerg.SUMMER_April-Sept.CONSISTENCY.precip_composite.daily_rainfall_bias_from_daily_climatology.BSISO_phase_"+str(BSISO)+"."+str(years[0])+"-"+str(years[-1])+"."+str(no_days).zfill(4)+"-days.lowerres_TRMM_grid.txt",  consistency[:,:], '%.4f')
	
	
			
	
			
			
		



