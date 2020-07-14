from netCDF4 import Dataset
import numpy as np
import glob
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.cm as cm
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, Normalize #DivergingNorm
import matplotlib.colors as mpcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import scipy as sp

def round_up_to_even(f):
	return math.ceil(f / 2.) * 2
	
#class MidpointNormalize(matplotlib.colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	#def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		#self.midpoint = midpoint
		#matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

	#def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		#x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		#return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
		
		
class MidpointNormalize(mpcolors.Normalize):
    def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False):
        self.vcenter = vcenter
        mpcolors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))
		


def map_composite_data(data, lats, lons, outfile, region, cmap, bounds, label,contour_or_grid,dtype=None):
	"""Plots a map of data (e.g. track density, precip) in the Southern Hemisphere"""
		
	if region == "WESTPAC":
		lat1=30
		lat2=0
		lon1=90
		lon2=140

	fig = plt.figure(figsize=(6, 3))
	ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])


	m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')

	m.drawcoastlines(linewidth=0.6, color='k') #gray
	m.drawcountries(linewidth=0.6, color='k') #gray

	cmap = plt.get_cmap(cmap)
	colors = cmap(np.linspace(0.1, 1.0, cmap.N))
	cmap2 = LinearSegmentedColormap.from_list(cmap, colors)

	#max = round_up_to_even(np.max(data))
	#print max
	#bounds = np.linspace(0, max, 11)
	#bounds = np.linspace(0,4000,9)
	print bounds
	#bounds = np.linspace(0, np.max(data))
	#norm = BoundaryNorm(bounds, ncolors=cmap2.N, clip=True)
	norm = BoundaryNorm(bounds, ncolors=cmap2.N)

	# mask 0 values so that we can set them to white/transparent rather than pale blue
	data = np.ma.masked_where(data == 0, data)
	lons, lats = np.meshgrid(lons, lats)

	#im = m.pcolormesh(lons, lats, data, shading='flat', cmap=cmap, latlon=True, norm=norm)
	#im.cmap.set_bad('white', alpha=0)  # set missing values to white/transparent
	
	
	midnorm = MidpointNormalize(vmin=bounds[0],vcenter=5., vmax=bounds[-1])
	print bounds[0]
	print bounds[-1]

	if dtype == 'bias':
		#if contour_or_grid == "grid":
		
		im = m.pcolormesh(lons, lats, data, shading='flat', cmap=cmap, latlon=True, norm = midnorm ) #vmin=bounds[0], vmax=bounds[-1])

		im.cmap.set_bad('silver') 
		cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="max") #,ticks=bounds
		#cbar.ax.set_yticklabels([i for i in bounds])

	elif dtype == "perc":
		im = m.pcolormesh(lons, lats, data, shading='flat', cmap=cmap2, latlon=True, norm=norm)
		#im.cmap.set_bad('grey')  # set missing values to white/transparent
		im.cmap.set_under('silver')
		#max_color=cmap(1.0)
		#im.cmap.set_over(max_color)
		cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, ticks=bounds,extend="max") #, extend="max",ticks=bounds[0::2]
		#for i in range(len(bounds)):
			#bounds[i] = int(bounds[i])
		#cbar.ax.set_yticklabels(bounds) #bounds[0::2]

	elif dtype == None:
		if contour_or_grid == "grid":
			im = m.pcolormesh(lons, lats, data, shading='flat', cmap=cmap2, latlon=True, norm=norm)
		elif contour_or_grid == "contour":
			im = m.contourf(lons, lats, data, cmap=cmap, latlon=True,levels=bounds,extend="max")
		im.cmap.set_bad('white', alpha=0)  # set missing values to white/transparent
		cbar = fig.colorbar(im, ax=ax,fraction=0.046, pad=0.04,extend="max")  # auto colorbar  #, extend="max"

	cbar.ax.tick_params(labelsize=12)  # colorbar font size
	plt.text(0.01, 1.02, label, transform=ax.transAxes, fontsize=10)
	plt.savefig(outfile, bbox_inches='tight', pad_inches=0.1, dpi=400)
	plt.close()



imerg_file_for_lon_lats = "/gws/nopw/j04/klingaman/emerton/GPM-IMERG/GPM_IMERG_3B-DAY.MS.MRG.3IMERG.2001.DAILY.V06_regridded_to_TRMM_grid.nc"
ff_imerg = Dataset(imerg_file_for_lon_lats,'r')
		
imerg_lons = ff_imerg.variables['longitude'][:]
imerg_lats = ff_imerg.variables['latitude'][:]


	
	
	
for month, mlabel in zip([1,2,3,4,5,6,7,8,9,10,11,12],['January','February','March','April','May','June','July','August','September','October','November','December']):
	
	IMERG_CLIM_DIR = "/gws/nopw/j04/klingaman/emerton/GPM-IMERG/"
	
	IMERG_CLIM_FILE = IMERG_CLIM_DIR+"GPM_IMERG.DAILY.TRMM_GRID.2001-2019.m"+str(month).zfill(2)+".95th_percentile.nc"
#
	print IMERG_CLIM_FILE
#
	ffCLIM = Dataset(IMERG_CLIM_FILE,'r')
#
	CLIM = ffCLIM.variables['precipitationCal'][:]
#
	ffCLIM.close()
	
	print CLIM
	
	print np.shape(CLIM)
	
	data = np.transpose(CLIM[0,:,:])
	
	print np.shape(data)
	
	outfile = "gpm-imerg.95th_percentile_of_climatology_precip.2001-2019.TRMM_grid.m"+str(month).zfill(2)+".png"
	
	label = "95th Percentile of Precipitation -  "+mlabel 
	bounds = np.linspace(0, 100, 11)
	map_composite_data(data, imerg_lats, imerg_lons, outfile, "WESTPAC", 'Blues', bounds, label,'grid',dtype="perc")
	
	print outfile
	

