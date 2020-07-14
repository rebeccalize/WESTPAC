# WESTPAC

This repository contains code used for the WESTPAC project to:

- use cdo to calculate precipitation climatologies (daily rainfall and 95th percentile extreme rainfall) using GPM IMERG satellite rainfall data 
- compute the bias of the daily rainfall from climatology under different conditions (ENSO, MJO, BSISO phases and combinations of these)
- compute the probability of extreme rainfall under different conditions (ENSO, MJO, BSISO phases and combinations of these)
- plot maps of the daily rainfall bias and probabilities of extreme rainfall over SE Asia 

The majority of the code is written in python and should be fairly well commented; the computations of bias and probabilities output text files on a regular grid (global, using the TRMM grid). The climatologies are netcdf files. 

The 'maps' directory contains all of the maps produced as part of this analysis, indicating rainfall biases and extreme rainfall probabilities in different phases of ENSO, MJO and BSISO, and combinations of these, across SE Asia:

- 'bias' maps indicate the mean difference between the daily rainfall in a given ENSO/MJO/BSISO phase, and the climatological daily rainfall
- 'consistency' maps indicate the number of days that have the same sign as the mean bias (i.e. if the mean bias at a certain grid point is positive, how many days in the dataset also had a positive bias rather than the opposite result?)
- 'probability' maps indicate the percentage of days in a given ENSO/MJO/BSISO phase, on which the daily rainfall exceeded the 95th percentile of climatology (typically using monthly 95th percentile values, calculated across 2001 - 2019; some use monthly percentiles calculated across only El Nino or La Nina years, or other variations - this should be evident in the filenames of the pngs)
- the 95th percentile rainfall values for each month are also mapped

For ENSO, NOAA's publicly available ONI index is used to define years in which the 3-month running mean index exceeded +/- 1.0. If the index exceeded +/- 1.0 at any point, then the El Nino / La Nina event is included. There are 3 El Nino events (2002-2003, 2009-2010, 2015-2016) and 4 La Nina events (2007-2008, 2010-2011, 2011-2012, 2017-2018) in the dataset (2001-2019). Since ENSO events typically peak in winter, a 'year' is taken as August of the ENSO event development year, through to July of the following year. For example, there was an El Nino event from SON 2009 through to MJJ 2010, so all days from August 2009 to July 2010 are used. 

For the MJO and BSISO, only days on which the amplitude exceeded 1.0 are included. 
