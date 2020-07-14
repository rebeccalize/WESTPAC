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
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, Normalize
import matplotlib.colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl

def round_up_to_even(f):
	return math.ceil(f / 2.) * 2
	
class MidpointNormalize(matplotlib.colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

	#def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		#x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		#return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

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
	colors = cmap(np.linspace(0.1, 1, cmap.N))
	cmap2 = LinearSegmentedColormap.from_list(cmap, colors)

	#max = round_up_to_even(np.max(data))
	#print max
	#bounds = np.linspace(0, max, 11)
	#bounds = np.linspace(0,4000,9)
	print bounds
	#bounds = np.linspace(0, np.max(data))
	norm = BoundaryNorm(bounds, ncolors=cmap2.N, clip=True)
	#norm = BoundaryNorm(bounds, ncolors=cmap2.N)

	# mask 0 values so that we can set them to white/transparent rather than pale blue
	data = np.ma.masked_where(data == 0, data)
	lons, lats = np.meshgrid(lons, lats)

	#im = m.pcolormesh(lons, lats, data, shading='flat', cmap=cmap, latlon=True, norm=norm)
	#im.cmap.set_bad('white', alpha=0)  # set missing values to white/transparent

	if dtype == 'bias':
		if contour_or_grid == "grid":
			im = m.pcolormesh(lons, lats, data, shading='flat', cmap=cmap, latlon=True, norm = MidpointNormalize(midpoint=0., vmin=bounds[0], vmax=bounds[-1]) ) #vmin=bounds[0], vmax=bounds[-1])
		elif contour_or_grid == "contour":
			plt.gca().patch.set_color('lightgrey') 
			im = m.contourf(lons, lats, data, cmap=cmap, latlon=True, levels=bounds, extend="both")
			#if mask_yes_no == 'yes':
			#im = m.contourf(lons,lats,mask,levels=[0,10],colors='lightgrey')
		im.cmap.set_bad('white')  # set missing values to white/transparent - was lightgrey , can't remember why from picsea
		cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="both")
		#cbar.ax.set_yticklabels([i for i in bounds])

	elif dtype == "perc":
		im = m.pcolormesh(lons, lats, data, shading='flat', cmap=cmap2, latlon=True, norm=norm)
		im.cmap.set_bad('white', alpha=0)  # set missing values to white/transparent
		cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04) #, extend="max"
		for i in range(len(bounds)):
			bounds[i] = int(bounds[i])
		#cbar.ax.set_yticklabels(bounds[0::2])

	elif dtype == None:
		if contour_or_grid == "grid":
			im = m.pcolormesh(lons, lats, data, shading='flat', cmap=cmap2, latlon=True, norm=norm)
		elif contour_or_grid == "contour":
			im = m.contourf(lons, lats, data, cmap=cmap, latlon=True,levels=bounds,extend="max")
		im.cmap.set_bad('white', alpha=0)  # set missing values to white/transparent
		cbar = fig.colorbar(im, ax=ax,fraction=0.046, pad=0.04,extend="max")  # auto colorbar  #, extend="max"

	cbar.ax.tick_params(labelsize=12)  # colorbar font size
	plt.text(0.76, 0.85, label, transform=ax.transAxes, fontsize=10)
	plt.savefig(outfile, bbox_inches='tight', pad_inches=0.1, dpi=400)
	plt.close()



imerg_file_for_lon_lats = "/gws/nopw/j04/klingaman/emerton/GPM_IMERG_3B-DAY.MS.MRG.3IMERG.2001.DAILY.V06_regridded_to_TRMM_grid.nc"
ff_imerg = Dataset(imerg_file_for_lon_lats,'r')
		
imerg_lons = ff_imerg.variables['longitude'][:]
imerg_lats = ff_imerg.variables['latitude'][:]


	
	
	
for BSISO in [1,2,3,4,5,6,7,8]: #7

	#infile = "gpm-imerg.precip_composite.daily_rainfall_bias_from_daily_climatology.MJO_phase_"+str(MJO)+".2001-2019*grid.txt"
	
	infile = "gpm-imerg.SUMMER_April-Sept.precip_composite.daily_rainfall_bias_from_daily_climatology.BSISO_phase_"+str(BSISO)+".2001-2019*.txt"
	
	
	infilename = glob.glob(infile)[0] #this gets the filename with the wildcard, since we don't know the number of days in the filename as it varies
	
	pcp_bias = np.genfromtxt(infilename, dtype=float)
	pcp_bias = np.transpose(pcp_bias) #lons, lats the wrong way round
	
	print pcp_bias
	print np.nanmax(pcp_bias)
	print np.nanmean(pcp_bias)
	print np.nanmin(pcp_bias)
	
	
	#no_days = infilename[92:96] #get the number of days out of the filename string
	
	no_days=infilename[112:116]
	
	print no_days
	
	outfile = "gpm-imerg.SUMMER_April-Sept.mean_daily_precip_bias_from_daily_climatology_composite.BSISO_phase_"+str(BSISO)+".2001-2019.mmperday.TRMM_grid.png"
	
	label = "BSISO Phase "+str(BSISO)+"\nApril - Sept\n       "+str(no_days)  #""TC-related precipitation (mm) " + str(y1) + "-" + str(y2)
	bounds = np.linspace(-12, 12, 13)
	map_composite_data(pcp_bias, imerg_lats, imerg_lons, outfile, "WESTPAC", 'RdBu', bounds, label,'grid',dtype="bias")
	
	print outfile
	
	
	#infile_consistency = "gpm-imerg.CONSISTENCY.precip_composite.daily_rainfall_bias_from_daily_climatology.MJO_phase_"+str(MJO)+".2001-2019*grid_v2.txt"
	
	infile_consistency = "gpm-imerg.SUMMER_April-Sept.CONSISTENCY.precip_composite.daily_rainfall_bias_from_daily_climatology.BSISO_phase_"+str(BSISO)+".2001-2019*.txt"
	
	print infile_consistency
	
	infile_consistency_name = glob.glob(infile_consistency)[0]
	
	consistency = np.genfromtxt(infile_consistency_name, dtype=float)
	consistency = np.transpose(consistency)
	
	print consistency
	print np.nanmean(consistency)
	print np.nanmax(consistency)
	print np.nanmin(consistency)
	
	
	outfile_consistency = "gpm-imerg.SUMMER_April-Sept.mean_daily_precip_bias_from_daily_climatology_composite.CONSISTENCY.BSISO_phase_"+str(BSISO)+".2001-2019.mmperday.TRMM_grid.png"
	bounds_consistency = np.linspace(0,100,21)
	#map_composite_data(consistency, imerg_lats, imerg_lons, outfile_consistency, "WESTPAC", 'YlGn', bounds_consistency, label, 'grid', dtype="perc")
	

