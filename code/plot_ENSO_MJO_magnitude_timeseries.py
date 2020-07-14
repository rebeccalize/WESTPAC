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
		
print ONI_EN
print ONI_LN
print ONI_n




fig = plt.figure()
fig.set_size_inches(25,4)
ax = fig.add_subplot(111, label="1")
ax2 = fig.add_subplot(111,label="2",frame_on=False)
	
xlen=[]

colours=['#3D5E96'] #,'#30C5D2' '#33ABC3','#3791B4','#3A77A5','#3D5E96','#404487','#442A78','#471069'
	
for MJO,c in zip([5], colours): #1,2,3,4,5,6,7,8
	MJO_file = "/gws/nopw/j04/klingaman/datasets/MJO_INDICES/MJO_phase"+str(MJO)+".jan-dec_dmeans_ts.1979-2019.nc"
	ffMJO = Dataset(MJO_file,'r')
	MJOamp = ffMJO.variables['phase_ts'][:]
	MJOdatesnc = ffMJO.variables['time'][:]

	t_unit = ffMJO.variables['time'].units
	t_cal = ffMJO.variables['time'].calendar
	tvalue = num2date(MJOdatesnc,units=t_unit, calendar=t_cal)
	
	MJOdates = [i.strftime("%Y-%m-%d") for i in tvalue]
	
	plotdates = matplotlib.dates.datestr2num(MJOdates)
		
	
	ax.bar(plotdates, MJOamp, width=1, color=c)
	ax.xaxis_date()
	


ax2.plot(range(len(ONI_EN)),ONI_EN,color='r')
ax2.plot(range(len(ONI_LN)),np.absolute(ONI_LN),color='b')
ax2.plot(range(len(ONI_n)),np.absolute(ONI_n),color='k')


ax2.yaxis.tick_right()
ax2.set_ylabel('abs(ONI)', color='blue',fontsize=14)
ax2.yaxis.set_label_position('right')
ax2.tick_params(axis='y',colors='red')
ax2.set_xticks([])

ax.set_xlim(np.min(plotdates),np.max(plotdates))
ax2.set_xlim(0,len(ONI_EN))
			
#print "xlen: ", xlen
#plt.xticks(fontsize=12)
#plt.xlim(0,65)
#xticklocs=np.arange(0,65, step=4)
#plt.ylim(0,300)


	
	
#plt.xticks(xticklocs, xticklocs/4) #sets the location of the xticks (one tick per day = every 4 timesteps), and the values (want to give it in days not timesteps, so /4)
#plt.yticks(fontsize=12)
#plt.xlabel('Time', fontsize = 14)
ax.set_ylabel('MJO Amplitude', color='#30C5D2',fontsize=14)
ax.tick_params(axis='y', colors='#30C5D2')

plt.savefig("ENSO_MJO_phases_5_amplitude_timeseries_1979-2020_test.png", bbox_inches='tight', pad_inches=0.05, dpi=500)
plt.close()



