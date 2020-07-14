#!/bin/bash

#This is not an extensive script computing all the different climatologies,
#I would just use it to run a set of cdo commands and then modify it to compute something slightly different when needed
#But this particular script shows the typical commands used to compute the monthly percentile values, in this case across only the El Nino or La Nina events
#Other commands used include cdo's 'ydaypctl' which calculated daily 95th percentile values across all the years, or 
#cdo's 'timpctl' command (as below), but across all days of all years (just one 95th percentile value per grid point), but 
#these were not used in the final analysis that resulted from various discussion meetings, and the method used should be evident in the filenames

#BSUB -W 12:00

cd GPM-IMERG
	
for m in 1 2 3 4 5 6 7 8 ; do

	cdo mergetime GPM_IMERG_DAILY.TRMM_grid.20030${m}.nc GPM_IMERG_DAILY.TRMM_grid.20100${m}.nc GPM_IMERG_DAILY.TRMM_grid.20160${m}.nc GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2003-2010-2016.m0${m}.nc
	
	cdo mergetime GPM_IMERG_DAILY.TRMM_grid.20080${m}.nc GPM_IMERG_DAILY.TRMM_grid.20110${m}.nc GPM_IMERG_DAILY.TRMM_grid.20120${m}.nc GPM_IMERG_DAILY.TRMM_grid.20180${m}.nc GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2008-2011-2012-2018.m0${m}.nc
	
	cdo timpctl,95 GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2003-2010-2016.m0${m}.nc -timmin GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2003-2010-2016.m0${m}.nc -timmax GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2003-2010-2016.m0${m}.nc GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2003-2010-2016.m0${m}.95th_percentile.nc
	
	cdo timpctl,95 GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2008-2011-2012-2018.m0${m}.nc -timmin GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2008-2011-2012-2018.m0${m}.nc -timmax GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2008-2011-2012-2018.m0${m}.nc GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2008-2011-2012-2018.m0${m}.95th_percentile.nc
	
done
	
for m in 9 ; do

	cdo mergetime GPM_IMERG_DAILY.TRMM_grid.20020${m}.nc GPM_IMERG_DAILY.TRMM_grid.20090${m}.nc GPM_IMERG_DAILY.TRMM_grid.20150${m}.nc GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2002-2009-2015.m0${m}.nc
	
	cdo mergetime GPM_IMERG_DAILY.TRMM_grid.20070${m}.nc GPM_IMERG_DAILY.TRMM_grid.20100${m}.nc GPM_IMERG_DAILY.TRMM_grid.20110${m}.nc GPM_IMERG_DAILY.TRMM_grid.20170${m}.nc GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2007-2010-2011-2017.m0${m}.nc
	
	cdo timpctl,95 GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2002-2009-2015.m0${m}.nc -timmin GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2002-2009-2015.m0${m}.nc -timmax GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2002-2009-2015.m0${m}.nc GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2002-2009-2015.m0${m}.95th_percentile.nc
	
	cdo timpctl,95 GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2007-2010-2011-2017.m0${m}.nc -timmin GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2007-2010-2011-2017.m0${m}.nc -timmax GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2007-2010-2011-2017.m0${m}.nc GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2007-2010-2011-2017.m0${m}.95th_percentile.nc
	
done
	
	
for m in 10 11 12 ; do

	cdo mergetime GPM_IMERG_DAILY.TRMM_grid.2002${m}.nc GPM_IMERG_DAILY.TRMM_grid.2009${m}.nc GPM_IMERG_DAILY.TRMM_grid.2015${m}.nc GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2002-2009-2015.m${m}.nc
	
	cdo mergetime GPM_IMERG_DAILY.TRMM_grid.2007${m}.nc GPM_IMERG_DAILY.TRMM_grid.2010${m}.nc GPM_IMERG_DAILY.TRMM_grid.2011${m}.nc GPM_IMERG_DAILY.TRMM_grid.2017${m}.nc GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2007-2010-2011-2017.m${m}.nc
	
	cdo timpctl,95 GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2002-2009-2015.m${m}.nc -timmin GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2002-2009-2015.m${m}.nc -timmax GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2002-2009-2015.m${m}.nc GPM_IMERG_DAILY.TRMM_grid.ElNinoYears.2002-2009-2015.m${m}.95th_percentile.nc
	
	cdo timpctl,95 GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2007-2010-2011-2017.m${m}.nc -timmin GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2007-2010-2011-2017.m${m}.nc -timmax GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2007-2010-2011-2017.m${m}.nc GPM_IMERG_DAILY.TRMM_grid.LaNinaYears.2007-2010-2011-2017.m${m}.95th_percentile.nc
	
done
