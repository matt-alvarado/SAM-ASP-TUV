#!/usr/bin/env python

# Based on build_ncep.py written by Robin Stevens, Dec. 2011
# Adapted by Chantelle Lonsdale, June 2017

# This script will download the necessary files for building a case using
# NCEP reanalysis-2 data or NARR data. Downloading it all can take a long
# time (perhaps hours), please be patient. They are very big files.

# This script will then go on to write the sfc, lsf, and snd files. 
# Reading the NCEP data takes time (perhaps 10 min?), so expect to wait for
# this as well. Keeping verbose mode off in the parametres below will
# shorten this time somewhat.

# If the location for the data you need is not inside the domain of the NARR 
# data, you will need to set NA to 0 below.

# Revision History
# edited by A. Dayalu 9/23/2021 to streamline new SAM-ASP_tuv run 
# (ie. new time and/or location)
#---------------------------------------------------------------------

import numpy as np
from datetime import *
import os
from time import localtime, strftime
import netCDF4

# --YOU WILL LIKELY ONLY NEED TO EDIT THIS SECTION---# 
# Set your input conditions in the params.txt file in the rundir that you
# set below!
# Default: the directory holding this script (the SAM-ASP-TUV checkout)
rundir = os.path.dirname(os.path.abspath(__file__)) #<==EDIT THIS if needed
caseidfill = "ADTEST" #<==EDIT THIS
#---- END EDIT SECTION ---#

inputfilename = "params.txt" #AD Test
# --- parameters you may need to change in params.txt file-----------------
#New instructions AD 9/23/2021
#In theory you shouldn't have to do anything here beyond changing 
#rundir above. 
#All other inputs would be in the params.txt file in your rundir.
#---------------------------------------------------------------------
#update the params.txt file (or create one) in the run directory (rundir)
# of your SAM-ASP run:

#The ','-sep param.txt file should look exactly like this (without the comment hashes)
#starting with the header, and case info on a new line 
# ** Note, in theory additional cases can be added to new lines, but this does not 
# seem like it is a currently available feature!!**
#lat,lon,year,month,day,hour,minute,Massemis,
#39.69583,-122.6,2009,12,17,10,0,4.86e-06,
# The only entries you will likely need to change are the
# latitude (degrees) [30,70]
# longitude (degrees) [-125,-75]
# Y,M,D,H.

# Note instructions later in this code for runs outside North America.

#--**old version of params.txt** had entries that were not used:--#
#lat,lon,time,Dpemis,Massemis,sigmaemis,zlower,zwidth,distance,xextent
#39.69583,-122.6,351.42,70.69247775,4.86e-06,1.14958271,1.13950504,0.89222415,5.64502954,1.54594687
#--end old version of params.txt--#

#--OLD instructions, out of date. Kept for backup. Don't Use:--
# text file containing tab-delimited parameters for the model runs.
# parameters should be in the following order:
# 3 condensation sink multiplier [0,5]
# 4 downward shortwave radiative flux (W/m2) [0,1165]
# 5 latitude (degrees) [30,70]
# 6 longitude (degrees) [-125,-75]
# 7 time coordinate (3-hourly data point from NCEP NARR. 0-247 are
#   points in January, 248-495 are points in July) [0-495]
#--END old instructions --#

#aerfilename = '12-10-04_GC_dNdlogD'

# a place to save the parameters for the model runs, with boundary
# layer heights and starting days added.
#extinputfilename = '13-02-22_GC_extinputs.txt'

# caseid to identify the data from this set of runs.
#caseidfill = 'W1101b'
#caseidfill = "ADTEST"

# directroy of your SAM-ASP run:
#rundir = '/rotor/data/p2033/SAM-ASP_ozone/'
#rundir = '/nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv'

# -- parameters for reading NCEP data --

# I'm going to write this assuming that we are inside NARR domain, but I'll leave in the options I had written for the build_ncep.py script anyway, in case we need them later.
NA = 1 # inside North America (domain of NARR)?
# 0 for false (use NCEP reanalysis-2 data)
# 1 for true (use NARR data)

p0 = 1000 # reference pressure (mb) for calculating potential temperature

npts = 3 # number of times to save to the lsf and snd files
# NOTE DIFFERENT FUNCTIONALITY FOR THIS SCRIPT FROM build_ncep.py:
# We want constant meteorology for the parameterization, so conditions
# will not evolve with time. SAM still has problems if files are not
# written such that model time does not exceed the range of time in
# the inputs, so we need more than one point.
# NCEP NARR data points are each 3 hours apart, so the last sounding written
# will be (npts-1)*3 hours after the startdate, if you are using NARR data.
# NCEP data points are each 6 hours apart, so the last sounding written
# will be (npts-1)*6 hours after the startdate for NCEP data.

samx = 'Ay' # switch for rotating winds.
# In NCEP data, uwind is zonal wind, vwind is meridional wind.
# In SAM, x should generally be the direction of the predominant wind.
# 'U': unidirectional wind. The wind flows in the model x direction.
#     (SAM uwind)=square root((NCEP uwind)**2 + (NCEP vwind)**2), (SAM vwind)=0
# 'N': x is North. (SAM uwind) = (NCEP vwind), (SAM vwind) = -(NCEP uwind)
# 'S': x is South. (SAM uwind) = -(NCEP vwind), (SAM vwind) = (NCEP uwind)
# 'E': x is East. (SAM uwind) = (NCEP uwind), (SAM vwind) = (NCEP vwind)
# 'W': x is West. (SAM uwind) = -(NCEP uwind), (SAM vwind) = -(NCEP vwind)
# 'T': x is theta radians (counter-clockwise) from East.
#     (SAM uwind) =  (NCEP uwind)*cos(theta) + (NCEP vwind)*sin(theta)
#     (SAM uwind) = -(NCEP uwind)*sin(theta) + (NCEP vwind)*cos(theta)
# 'Ax': x is autofitted to the mean wind direction in the >900 hPa layer
#     (roughly the boundary layer) for each 'sounding'. Note that this
#     means that x could be North at one time and North-West later.
# 'Ay': y is autofitted to the mean wind direction in the >900 hPa layer
#     for each 'sounding'. Useful for 2D wall model.
# If samx un
# If samx undefined, use 'U' setting.

theta = np.pi*0.40
verbose = 0

# print extra output for checking that latitude, longitude, time is
# consistent and correct for all files read. 0 = false, 1 = true.

cl = '-nc' # overwrite files?
# set to '-nc' to not download files if one with the same name already exists.
# set to '-r' to download anyway and overwrite. (I think this works, not sure)

# NARR data for North America
narrurl = 'ftp://ftp.cdc.noaa.gov/Datasets/NARR/'
# NCEP reanalysis data - global, but much lower resolution
if NA:
   ncepurl = 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/'
else:
   ncepurl = 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis2/'

# -- parameters for queueing SAM: ---
samExecName = './SAM_RAD_CAM_MICRO_SAM1MOM'

maxn = 25000 # maximum number of timesteps for model to run

# if maxd < dd, then (I think) only one output will be generated for each case,
# at the timestep corresponding to ddes.
maxd = 10000. # maximum distance from source for which data is wanted [m]
mind =  5000. # minimum distance from source for which data is wanted [m]
dd   = 50000. # difference in distance between successive model output. [m]

# Note that the minimum windspeed for which maxd will be achievable is
#  maxd/(maxn*dt)
# For lower windspeeds, the model will stop short of this.
# The minimum windspeed for which mind will be achievable is
#   mind/(maxn*dt)
# For lower windspeeds, no model run with be started.

dxy = '500.'
noxshift = '1.0'
temis2D = 1800.

dt = 8. # model timestep [s]

mode = 'CEM' # Cloud Resolving Model mode. LES not valid for 100+m boxes
wetness = '0.4' # 0 = absolutely dry surface, 1 = water surface

nuc = ('0','0','0') # nucleation settings as (bin_nuc, tern_nuc, ion_nuc) 
ionrate = '0.0' # ionization rate, only used with ion_nuc

tauls = '3600' # nudging timescale

#emisgridx = '64'
#emisgridy = '1'

# background gases (and things that control them)
#nh3init = '0.0' # intended background NH3 [ppt]
#nh4so4rat = '0.0' # intended NH4:SO4 ratio

# -- constants --

rc = 8.314472/29.19 # ratio of gas constant to specific heat capacity of air

airdens = 1.2 # density of air, in kg/m^3

HoverT = 8.31432/(0.0289644*9.80665) # scale height divided by temperature [m/K]

# --- end parameters -------------------------------------------------


def fillarrays(infile,varstring):
   lat  =infile.variables["lat"][:]   # latitude
   lon  =infile.variables["lon"][:]   # longitude
   time =infile.variables["time"][:]  # time (hours since 1800-1-1 or 0001-1-1)

   #offset = infile.variables[varstring].add_offset # offset for data
   #scalef = infile.variables[varstring].scale_factor # scaling factor for data

   return([lat,lon,time])

def findpoint(lat,lon,time,emislat,emislon,ncephours1800,ncephours0001,NA):
   # find the x,y,t values that correspond to the lat,lon,time values we
   # are looking for.

   # find the x and y values for the lat and lon we want
   bestx = 0
   besty = 0
   bestfit = 9999.9
   if NA:
      for x in range(0,349):
         for y in range(0,277):
            thisfit = (lat[y,x]-emislat)**2 + (lon[y,x]-emislon)**2
            if thisfit < bestfit:
               bestfit = thisfit
               bestx = x
               besty = y
   else:
      for x in range(0,len(lon)):
         for y in range(0,len(lat)):
            thisfit = (lat[y]-emislat)**2 + (lon[x]-emislon)**2
            if thisfit < bestfit:
               bestfit = thisfit
               bestx = x
               besty = y

   # find the frame that contains the date and time we want
   bestfit = 9999.9
   bestt = 0
   if np.median(time[:]) < 1E7:
      for t in range(np.shape(time)[0]):
         thisfit = abs(time[t]-ncephours1800)
         if thisfit < bestfit:
            bestfit = thisfit
            bestt = t
      t = bestt
   else:
      for t in range(np.shape(time)[0]):
         thisfit = abs(time[t]-ncephours0001)
         if thisfit < bestfit:
            bestfit = thisfit
            bestt = t
      t = bestt

   return([bestx,besty,t])

def printstuff(x,y,t):
   print 'i = '+str(x)+' j = '+str(y)
   print str(emisdate)
   print 't:'+str(t)
   print 'time used:'+str(time[t])
   return

def printstuff2(x,y,t):
   if NA:
      print 'lat = ' + str(infile.variables["lat"][y,x])
      print 'lon = ' + str(infile.variables["lon"][y,x])
   else:
      print 'lat = ' + str(infile.variables["lat"][y])
      print 'lon = ' + str(infile.variables["lon"][x])
   print str(emisdate)
   print 't:'+str(t)
   print 'time used'+str(time[t])
   return

def number_dist(Dpm,sigma,mass_flux,aerodens):
   '''Change a mass flux, kg/m2/s to #/second'''
   mass_box=mass_flux*float(dxy)**2 #kg/s
   N0=1000. # this is just some dummy number for now
   (xk,x,dpk,dp) = TOMAS_arrays()
   dist_new=logdist(N0,Dpm*1E-9, sigma, dp*1E-6) # put this in tomas.py file from c.py, check units here
   mass_array=4/3.*((dp*1E-6)/2.)**3*aerodens #volume of particle in bin times density
   dist=sizedist(dist_new,'logDp')
   dist.conc(dpk*1E-6)
   total_mass=(mass_array*dist.con).sum()  
   div=mass_box/total_mass
   numbers=N0*div # number total per box per secon
   dist_new=logdist(numbers,Dpm*1E-9,sigma,dp*1E-6)
   dist=sizedist(dist_new,'logDp')
   dist.conc(dpk*1E-6)
   total_mass=(mass_array*dist.con).sum()
   return numbers/(float(dxy)**2*dz*1E6) # number concentration [per cm-3] ,total_mass, mass_box, mass_array

# --------------------------------------------------------------------
# 1) Read the inputs file, to see which paremeters we are using for
# this model run.
# --------------------------------------------------------------------

inputfile = open(inputfilename,'r')
lines = inputfile.readlines()
inputfile.close()
header = lines.pop(0)

ncases = len(lines)

#aerfile = open(aerfilename,'r')
#aerlines = aerfile.readlines()
#aerfile.close()
#header = aerlines.pop(0)



#if len(aerlines) != ncases:
#   print('len(aerlines) != ncases')
#   print('len(aerlines)',len(aerlines),'ncases',ncases)
#   print('quitting...')
#   quit()

caseid = 1

for lineid in range(0,1):#ncases): # here the looping begins if you have multiple days!

   caseidstr = '%(cid)04d' % {"cid":caseid}
   print ' case '+caseidstr+'/%(len)3d' % {"len":len(lines)}

   
   # unpack the parameters
   vals=lines[lineid].split(',')
   vals[-1]=vals[-1].split()
   emislat = float(vals[0]) # latitude (degrees)
   emislon = float(vals[1]) # longitude (degrees)
   input_year = int(vals[2])
   input_month = int(vals[3])
   input_day = int(vals[4])
   input_hour = int(vals[5])
   input_minute = int(vals[6])
   #calculate decimal day of year (day 1 = Jan 1 of year); this will be
   #used later and shuttled to variable timepoint
   timepoint = (datetime(input_year, input_month, input_day, input_hour, input_minute)-datetime(input_year, 1, 1)).total_seconds()/86400+1

   #timepoint = int(round(float(vals[2]))) # time coordinate
   massemis = float(vals[7])

   #Not used for now ....
   #xextent = int(round(float(vals[8])*1000./float(dxy))) # form km to m to # of boxes
   #yextent = xextent
   #zlower = int(round(float(vals[6])))
   #zextent = int(round(float(vals[7])))
   #zhigher = zlower+zextent

   #casename = ('SAMASP'+caseidstr)
   casename = (caseidfill)
   os.system('mkdir -p RUNS/'+casename)

# --------------------------------------------------------------------
# 2) Set up the case for this model run. This means writing the sfc,
# lsf, and snd files, as well as the prm file.
# --------------------------------------------------------------------

   startdate = datetime(input_year,input_month,input_day,input_hour,input_minute)

   #timepoint
   print('startdate = ',startdate)   
   if NA == 0 and emislon < 0.:
      emislon = emislon+360.

   # silly strings!
   yrstr = str(startdate.year)
   if startdate.month < 10: mnstr = '0'+str(startdate.month)
   else: mnstr = str(startdate.month)

   nceptime1800 = startdate - datetime(1800,1,1) # timedelta object, time since 1800
   ncephours1800 = nceptime1800.days*24+nceptime1800.seconds/3600 # time since 01/01/1800 in hours
   nceptime0001 = startdate - datetime(1,1,1) # timedelta object, time since year 1
   ncephours0001 = nceptime0001.days*24+nceptime0001.seconds/3600 # time since 01/01/01 in hours

   # --- download the data we need --------------------------------------
   
   if NA: # get NARR data
      # get the air temperature file
      p = '-P '+rundir+'/NCEP/'
      os.system('wget '+cl+' '+ p +' '+narrurl+'pressure/air.'+yrstr+mnstr+'.nc')
      # get the specific humidity file
      os.system('wget '+cl+' '+ p +' '+narrurl+'pressure/shum.'+yrstr+mnstr+'.nc')
      # get the surface pressure file
      os.system('wget '+cl+' '+ p +' '+narrurl+'monolevel/pres.sfc.'+yrstr+'.nc')
      # get the planetary boundary layer height file
      os.system('wget '+cl+' '+ p +' '+narrurl+'monolevel/hpbl.'+yrstr+'.nc')
      # get the u-wind file
      os.system('wget '+cl+' '+ p +' '+narrurl+'pressure/uwnd.'+yrstr+mnstr+'.nc')
      # get the v-wind file
      os.system('wget '+cl+' '+ p +' '+narrurl+'pressure/vwnd.'+yrstr+mnstr+'.nc')
      # --- sfc files: ---
      # get the sensible heat flux file
      os.system('wget '+cl+' '+ p +' '+narrurl+'monolevel/shtfl.'+yrstr+'.nc')
      # get the latent heat flux file
      os.system('wget '+cl+' '+ p +' '+narrurl+'monolevel/lhtfl.'+yrstr+'.nc')
      # get the u-momentum flux file
      # no narr file for this2
      os.system('wget '+cl+' '+ p +' '+ncepurl+'surface_gauss/uflx.sfc.gauss.'+yrstr+'.nc')
      # get the v-momentum flux file
      # no narr file for this
      os.system('wget '+cl+' '+ p +' '+ncepurl+'surface_gauss/vflx.sfc.gauss.'+yrstr+'.nc')
   else: # get NCEP reanalysis-2 data
      os.system('wget '+cl+' '+ p +' '+ncepurl+'pressure/air.'+yrstr+'.nc')
      # No NCEP reanalysis-2 4/day data for specific humidity
      os.system('wget '+cl+' '+ p +' '+'ftp://ftp.cdc.noaa.gov/Datasets/'+\
                   'ncep.reanalysis/pressure/shum.'+yrstr+'.nc')
      os.system('wget '+cl+' '+ p +' '+ncepurl+'surface/pres.sfc.'+yrstr+'.nc')
      os.system('wget '+cl+' '+ p +' '+ncepurl+'pressure/uwnd.'+yrstr+'.nc')
      os.system('wget '+cl+' '+ p +' '+ncepurl+'pressure/vwnd.'+yrstr+'.nc')
      os.system('wget '+cl+' '+ p +' '+ncepurl+'gaussian_grid/shtfl.sfc.gauss.'+yrstr+'.nc')
      os.system('wget '+cl+' '+ p +' '+ncepurl+'gaussian_grid/lhtfl.sfc.gauss.'+yrstr+'.nc')
      os.system('wget '+cl+' '+ p +' '+ncepurl+'gaussian_grid/uflx.sfc.gauss.'+yrstr+'.nc')
      os.system('wget '+cl+' '+ p +' '+ncepurl+'gaussian_grid/vflx.sfc.gauss.'+yrstr+'.nc')
   
   # --- read the data files for writing the lsf and snd files ----------

   # --- get air temperature data ---
   print 'getting air temperature data...'
   if NA:
      infile = netCDF4.Dataset('NCEP/air.'+yrstr+mnstr+'.nc','r',format='netcdf3') # Read file using NIO
   else:
      infile = netCDF4.Dataset('NCEP/air.'+yrstr+'.nc','r',format='netcdf3') # Changed for NCEP date coding

   level=infile.variables["level"][:] # pressure level (mb)

   [lat,lon,time] = fillarrays(infile,"air")

   [bestx,besty,t] = findpoint(lat,lon,time,emislat,emislon,ncephours1800,ncephours0001,NA)

   # find day of year, including fraction
   # For this queueing script, meteorological variables are not
   # changing with time. Therefore, the number of hours between points
   # is irrelevent. Let's choose 1 day, so that we need not print many
   # data points to the files.
   if NA:
      emisdate = datetime(1800,1,1,0) + timedelta(hours=time[t])
   else:
      emisdate = datetime(1,1,1,0) + timedelta(hours=time[t])
   day = np.zeros(npts)
   ofyear = emisdate - datetime(emisdate.year,01,01)
   case_year = emisdate.year # calendar year of day0 (year0 in prm)
   for i in range(npts):
      day[i] = ofyear.days + float(ofyear.seconds)/(3600.*24.) + i


   # the 3D data from the output files comes out as data[time,level,y,x]
   # air temperature (K)
   temp = infile.variables["air"][t,:,besty,bestx]  
   #Above strings all changed from t-1 to t for lower bounds.
   infile.close()

   # convert temperature to potential temperture
   ptemp = temp*(p0/level[:])**rc


   # --- get the surface pressure ---
   print 'getting surface pressure...'
   infile = netCDF4.Dataset('NCEP/pres.sfc.'+yrstr+'.nc','r')

   [lat,lon,time] = fillarrays(infile,"pres")

   [bestx,besty,t] = findpoint(lat,lon,time,emislat,emislon,ncephours1800,ncephours0001,NA)
   
   pres = infile.variables["pres"][t,besty,bestx]  
   pres = pres*0.01 # convert from Pa to mb or hPa
   infile.close()

   # --- get planetary boundary layer height data ---
   print 'getting bl height data...'
   # I don't have a plan for this outside NA at the moment
   infile=netCDF4.Dataset('NCEP/hpbl.'+yrstr+'.nc','r',format='netcdf3')

   [lat,lon,time] = fillarrays(infile,"hpbl")

   [bestx,besty,t] = findpoint(lat,lon,time,emislat,emislon,ncephours1800,ncephours0001,NA)

   # height of the planetary boundary layer [m]
   hpbl = infile.variables["hpbl"][t,besty,bestx]  
   infile.close()

   # convert the hpbl into pressure space
   # I'm using the temperature just above the surface here, the resulting
   # errors should be less than 2%
   hpblp = pres*np.exp(-hpbl/(HoverT*temp[0]))

   # --- get u-wind data ---
   print 'getting u-wind data...'
   if NA:
      infile=netCDF4.Dataset('NCEP/uwnd.'+yrstr+mnstr+'.nc','r',format='netcdf3') 
   else:
      infile=netCDF4.Dataset('NCEP/uwnd.'+yrstr+'.nc','r',format='netcdf3') 


   [lat,lon,time] = fillarrays(infile,"uwnd")

   [bestx,besty,t] = findpoint(lat,lon,time,emislat,emislon,ncephours1800,ncephours0001,NA)

   # u-wind (m/s)
   uwind = infile.variables["uwnd"][t,:,besty,bestx]  
   infile.close()


   #--- get the v-wind data ---
   print 'getting v-wind data...'
   if NA:
      infile=netCDF4.Dataset('NCEP/vwnd.'+yrstr+mnstr+'.nc','r',format='netcdf3') 
   else:
      infile=netCDF4.Dataset('NCEP/vwnd.'+yrstr+'.nc','r',format='netcdf3') 

   [lat,lon,time] = fillarrays(infile,"vwnd")

   [bestx,besty,t] = findpoint(lat,lon,time,emislat,emislon,ncephours1800,ncephours0001,NA)

   # v-wind (m/s)
   vwind = infile.variables["vwnd"][t,:,besty,bestx]  
   infile.close()

   # rotate winds, if necessary
   if samx == 'N':
      uls = vwind
      vls = -uwind
   elif samx == 'S':
      uls = -vwind
      vls = uwind
   elif samx == 'E':
      uls = uwind
      vls = vwind
   elif samx == 'W':
      uls = -uwind
      vls = -vwind
   elif samx == 'T':
      uls =  uwind*np.cos(theta) + vwind*np.sin(theta)
      vls = -uwind*np.sin(theta) + vwind*np.cos(theta)
   elif samx == 'Ax': # autofit wind x
      hpblz = np.where(level > min(hpblp,pres-50))[0][0]
      wspd = (uwind[hpblz]**2+vwind[hpblz]**2)
      # Weight wind direction by wind speed
      theta = ( np.arctan2(vwind[hpblz],uwind[hpblz])*wspd/wspd.mean() ).mean()
      uls =  uwind*np.cos(theta) + vwind*np.sin(theta)
      vls = -uwind*np.sin(theta) + vwind*np.cos(theta)
   elif samx == 'Ay': # autofit wind y
      hpblz = np.where(level > min(hpblp,pres-50))[0][0]
      print(hpblp,pres,min(hpblp,pres-50),level)
      wspd = (uwind[hpblz]**2+vwind[hpblz]**2)
      # Weight wind direction by wind speed
      theta = ( np.arctan2(vwind[hpblz],uwind[hpblz])*wspd/wspd.mean() ).mean()
      uls = uwind*np.sin(theta) - vwind*np.cos(theta)
      vls = uwind*np.cos(theta) + vwind*np.sin(theta)
   else: # unidirectional wind. All wind flows in the model x direction.
      uls = np.sqrt(uwind**2+vwind**2)
      vls = vwind*0.0
   
   # find mean wind speed in the boundary layer
   vg = 0.0
   for z in range(len(level[:])):
      if level[z] > min(hpblp,pres-50.):
         vg = vg*z/(z+1.)+vls[z]/(z+1.)

   # write the extended inputfile
   #extinputfile = open(extinputfilename,'a')
   #extinputfile.write(lines[lineid][:-1]+'\t'+str(hpbl)+'\t'+str(day[0])+'\t'+str(vg)+'\n')
   #extinputfile.close()

   # --- get the specific humidity ---
   print 'getting specific humidity data...'
   if NA:
      infile=netCDF4.Dataset('NCEP/shum.'+yrstr+mnstr+'.nc','r',format='netcdf3')
   else:
      infile=netCDF4.Dataset('NCEP/shum.'+yrstr+'.nc','r',format='netcdf3')

   [lat,lon,time] = fillarrays(infile,"shum")

   [bestx,besty,t] = findpoint(lat,lon,time,emislat,emislon,ncephours1800,ncephours0001,NA)

   # specific humidity (kg/kg)
   shum = infile.variables["shum"][t,:,besty,bestx]  
   infile.close()

   # convert specific humidity to mixing ratio (g/kg)
   wmix = shum / (1.0 - shum) * 1000.0


   # --- write the snd file ---
   print 'writing snd file...'
   #outputfile =  open('/home/ksakamoto/SAM/RUNS/'+casename +'/snd','w')
   outputfile =  open(rundir+'/RUNS/'+casename+'/snd','w')
   outputfile.write(' z[m] p[mb] tp[K] q[g/kg] u[m/s] v[m/s]\n')
   for i in range(npts):
      outputfile.write(' '+str(day[i])+',  '+str(len(level[:]))+', '+str(pres)[:6]+'  day,levels,pres0\n')
      for z in range(len(wmix[:])):
         outputfile.write('-999.\t'+str(level[z])+'\t'+str(ptemp[z])[:6]+'\t'+str(wmix[z])[:6] +'\t'+str(uls[z])[:5]+'\t'+str(vls[z])[:5]+'\n')
   outputfile.close()


   # --- write the lsf file ---
   print 'writing lsf file...'
   outputfile =  open(rundir+'/RUNS/'+casename+'/lsf','w')
   outputfile.write(' z[m] p[mb]  tls[K/s] qls[kg/kg/s] uls  vls  wls[m/s]\n')
   for i in range(npts):
      outputfile.write(' '+str(day[i])+',  '+str(len(level[:]))+', '+str(pres)[:6]+'  day,levels,pres0\n')
      for z in range(len(wmix[:])):
         outputfile.write('-999.\t'+str(level[z])+'\t0.\t0.\t'+str(uls[z])[:5]+'\t'+str(vls[z])[:5]+'\t0.\n')
   outputfile.close()


   # --- read the data files for writing the sfc file -------------------

   # --- read the sensible heat flux ---
   print 'reading the sensible heat flux...'
   if NA:
      infile = netCDF4.Dataset('NCEP/shtfl.'+yrstr+'.nc','r',format='netcdf3')
   else:
      infile = netCDF4.Dataset('NCEP/shtfl.sfc.gauss.'+yrstr+'.nc','r',format='netcdf3') # Read file using NIO

   [lat,lon,time] = fillarrays(infile,"shtfl")

   [bestx,besty,t] = findpoint(lat,lon,time,emislat,emislon,ncephours1800,ncephours0001,NA)

   sens = infile.variables["shtfl"][t,besty,bestx]  
   infile.close()

   # --- read the latent heat flux ---
   print 'reading the latent heat flux...'
   if NA:
      infile = netCDF4.Dataset('NCEP/lhtfl.'+yrstr+'.nc','r',format='netcdf3')
   else:
      infile = netCDF4.Dataset('NCEP/lhtfl.sfc.gauss.'+yrstr+'.nc','r',format='netcdf3')


   [lat,lon,time] = fillarrays(infile,"lhtfl")

   [bestx,besty,t] = findpoint(lat,lon,time,emislat,emislon,ncephours1800,ncephours0001,NA)

   latent = infile.variables["lhtfl"][t,besty,bestx]  
   infile.close()


   # ----- read in the momentum fluxes -----
   
   if NA: # timestep different from other files, so we have to muck around

      # --- read the u-momentum flux ---
      print 'reading the u-momentum flux...'
      infile = netCDF4.Dataset('NCEP/uflx.sfc.gauss.'+yrstr+'.nc','r',format='netcdf3') 

      [lat,lon,time] = fillarrays(infile,"uflx")
      
      # find x and y that best represents lat and lon

      # Longitude here goes from 0 to 360, not -180 to 180
      if emislon < 0.:
         emislon1=emislon+360.
      else:
         emislon1=emislon

      bestx = 0
      while abs(lon[bestx+1]-(emislon1))<abs(lon[bestx]-(emislon1)):
         bestx+=1

      besty = 0
      while abs(lat[besty+1]-emislat)<abs(lat[besty]-emislat):
         besty+=1

      uflx = infile.variables["uflx"][:,besty,bestx]  
      infile.close()

      # --- read the v-momentum flux ---
      print 'reading the v-momentum flux...'
      infile = netCDF4.Dataset('NCEP/vflx.sfc.gauss.'+yrstr+'.nc','r',format='netcdf3') 

      [lat,lon,time] = fillarrays(infile,"vflx")

      vflx = infile.variables["vflx"][:,besty,bestx]  
      infile.close()

      # get scalar momentum flux, and calculate tau
      mflx_tmp = (uflx**2 + vflx**2)**0.5 / airdens

      # determine whether time used for mflx is the same as what we want
      mflx_1stdate = datetime(1,1,1,0) + timedelta(hours=time[0])
      mflx_1stofyear = mflx_1stdate - datetime(mflx_1stdate.year-1,12,31)
      mflx_1stday = mflx_1stofyear.days + float(mflx_1stofyear.seconds)/(3600.*24.)

      # find time[t] that closest to, but not less than day[i]
      t = 0
      mflx_day = mflx_1stday
      while mflx_day < day[0]:
         t+=1
         mflx_day = mflx_1stday+(time[t]-time[0])/24.
      # linearly interpolate between mflx points
      f = (mflx_day-day[0])/( (time[t]-time[t-1])/24. )
      mflx = mflx_tmp[t]*(1-f) + mflx_tmp[t-1]*f

   else: # we are using NCEP data, timestep same as other files

      # --- read the u-momentum flux ---
      print 'reading u-momentum flux...'
      infile = netCDF4.Dataset('NCEP/uflx.sfc.gauss.'+yrstr+'.nc','r',format='netcdf3') 

      [lat,lon,time] = fillarrays(infile,"uflx")

      # NCEP format, need to find lat, lon, time again.
      
      bestx = 0
      # Longitude here goes from 0 to 360, not -180 to 180
      if emislon < 0.:
         emislon1=emislon+360
      else:
         emislon1=emislon
      while abs(lon[bestx+1]-(emislon1))<abs(lon[bestx]-(emislon1)):
         bestx+=1

      besty = 0
      while abs(lat[besty+1]-emislat)<abs(lat[besty]-emislat):
         besty+=1

      bestfit = 9999.9
      bestt = 0
      if time[0] < 1E7:
         for t in range(0,len(time[:])):
            thisfit = abs(time[t]-ncephours1800)
            if thisfit < bestfit:
               bestfit = thisfit
               bestt = t
         t = bestt
         emisdate = datetime(1800,1,1,0) + timedelta(hours=time[t])
      else:
         for t in range(0,len(time[:])):
            thisfit = abs(time[t]-ncephours0001)
            if thisfit < bestfit:
               bestfit = thisfit
               bestt = t
         t = bestt
         emisdate = datetime(1,1,1,0) + timedelta(hours=time[t])

      if (verbose):
         print ' to check correct lat, lon, time in ulfx.sfc.gauss.'+yrstr+'.nc:'
         printstuff2(bestx,besty,t)

      uflx = infile.variables["uflx"][t,besty,bestx]  
      infile.close()

      # --- read the v-momentum flux ---
      print 'reading the v-momentum flux...'
      infile = netCDF4.Dataset('NCEP/vflx.sfc.gauss.'+yrstr+'.nc','r',format='netcdf3') 

      [lat,lon,time] = fillarrays(infile,"vflx")

      bestfit = 9999.9
      bestt = 0
      if time[0] < 1E7:
         for t in range(0,len(time[:])):
            thisfit = abs(time[t]-ncephours1800)
            if thisfit < bestfit:
               bestfit = thisfit
               bestt = t
         t = bestt
         emisdate = datetime(1800,1,1,0) + timedelta(hours=time[t])
      else:
         for t in range(0,len(time[:])):
            thisfit = abs(time[t]-ncephours0001)
            if thisfit < bestfit:
               bestfit = thisfit
               bestt = t
         t = bestt
         emisdate = datetime(1,1,1,0) + timedelta(hours=time[t])

      if (verbose):
         printstuff2(bestx,besty,t)

      vflx = infile.variables["vflx"][t,besty,bestx]  
      infile.close()

      # get scalar momentum flux, and calculate tau
      mflx = (uflx**2 + vflx**2)**0.5 / airdens
      

   # --- write the sfc file ---
   print 'writing sfc file...'
   outputfile = open(rundir+'/RUNS/'+casename+'/sfc','w')
   outputfile.write('day sst(K) H(W/m2) LE(W/m2) TAU(m2/s2)\n')
   for i in range(npts):
      outputfile.write(' '+str(day[i])+'\t0.\t'+str(-sens)[:6]+'\t'+str(-latent)[:6]+'\t'+str(mflx)[:6]+'\n')

   outputfile.close()


   # --- write the prm file ---
   print 'writing prm file...'

   # The casename and caseid are each limited to 40 characters.
   # I cannot fit all of the parameters in either one, so I am going
   # to identify the cases by a single number that refers to the line
   # of the inputs file that has the parameters for that case.
   # So long as the inputs file is unchanged, this will work fine.
   # If it is changed, there will still be clues in the runtime files.
   # Latitude, longitude, day of year, so2 emissions, and nox emissions
   # are all printed to the OUT_RUNTIME file.

   # This is the mass emis to number per box stuff for KS
   #Noemis = number_dist(Dpemis,sigmaemis,massemis,1400.)

   doSAM = True

   #ndes = int(desdist/(vg*dt))+int(temis2D/dt)+1
   #if ndes > maxn:
   #   print 'huge ndes, vg=',vg,'m/s, desdist=',desdist,'m'
   #   print 'skipping this run'
   #   doSAM = False

   nsave3Dstart = int(mind/(vg*dt))+int(temis2D/dt)+1
   #if nsave3Dstart > maxn:
   #   print 'huge nsave3Dstart, vg<'+str(mind/(maxn*dt))+' m/s'
   #   print 'skipping this run'
   #   print 'nsave3Dstart: '+str(nsave3Dstart)+' vg: '+str(vg)
   #   doSAM = False

   nstop = int(1800) # two hours #int(maxd/(vg*dt))+int(temis2D/dt)+1

   if doSAM:

      nsave3D = int(15*60/dt) # every 15 minutes int(dd/(vg*dt))
      # ensure ndes is a multiple of nsave3D
   #   if ndes < nsave3D: nsave3D = ndes
   #   else:
   #      while ndes%nsave3D>0: nsave3D+=1
      # ensure nstop is a multiple of nsave3D
      nstop = int(2*3600/dt) #int(1+nstop/nsave3D)*nsave3D

      if nstop > maxn:
         print 'huge nstop, vg<'+str(maxd/(maxn*dt))+' m/s'
         print 'setting nstop='+str(maxn)
         print 'nstop: '+str(nstop)+' vg: '+str(vg)
         nstop = maxn

      prmfile =  open(rundir+'/RUNS/'+casename +'/prm','w')
      prmfile.write(' &MICRO_M2005\n')
      prmfile.write(' doicemicro = .false.\n\n')  
             
      prmfile.write(' &PARAMETERS')
      prmfile.write(' caseid = \''+caseidfill+'\',\n')
      prmfile.write(' nrestart = 0,\n\n')

      prmfile.write(' '+mode+' = .true.,\n')
      prmfile.write(' LAND = .true.,\n')
      prmfile.write(' OCEAN = .false.,\n\n')
      #prmfile.write(' soil_wetness = '+wetness+'\n\n')

      prmfile.write(' dosgs        = .true.,\n')
      prmfile.write(' dodamping    = .true.,\n') # should be irrelevent
      prmfile.write(' doupperbound = .false.,\n') # should be irrelevent
      #prmfile.write(' dosmagor     = .false.,\n')
      #prmfile.write(' doscalar     = .false.,\n')
      prmfile.write(' docloud      = .true.,\n')
      prmfile.write(' doprecip     = .false., \n')
      prmfile.write(' dolongwave   = .false.,\n')
      prmfile.write(' doshortwave  = .false.,\n')
      prmfile.write(' dosurface    = .true.,\n')
      prmfile.write(' dolargescale = .true.,\n')
      prmfile.write(' doradforcing = .false.,\n')
      prmfile.write(' dotracers    = .true.,\n')
      prmfile.write(' dosfcforcing = .true.,\n')
      prmfile.write(' docoriolis   = .true.,\n')
      prmfile.write(' donudging_uv = .true.,\n')
      prmfile.write(' donudging_tq = .true.,\n')
      prmfile.write(' tauls        = '+tauls+'.,\n\n')

      prmfile.write(' SFC_FLX_FXD    = .true.,\n')
      prmfile.write(' SFC_TAU_FXD    = .true.,\n\n')

      prmfile.write(' latitude0 = '+str(emislat)+',\n')
      prmfile.write(' longitude0 = '+str(emislon)+',\n\n')

      prmfile.write(' fcor = 0.376e-4,\n\n')

      prmfile.write(' dx = '+dxy+',\n')
      prmfile.write(' dy = '+dxy+',\n')
      prmfile.write(' dt = '+str(dt)+',\n\n')

      prmfile.write(' day0= '+str(day[0])+',\n')
      # calendar year of day0 (TUV's date)
      prmfile.write(' year0 = '+str(case_year)+',\n\n')

      prmfile.write(' nstop = '+str(nstop)+',\n')
      prmfile.write(' nrestart_skip = '+str(nstop)+',\n\n')
      prmfile.write(' nprint = 100,\n')
      prmfile.write(' nstat = 100,\n')
      prmfile.write(' nstatfrq = 5,\n\n')

      prmfile.write(' doSAMconditionals = .false.\n')
      prmfile.write(' dosatupdnconditionals = .false.\n\n')

      prmfile.write(' nsave2D      = 20,\n')
      prmfile.write(' nsave2Dstart	= 99999999,\n')
      prmfile.write(' nsave2Dend   = 999164000,\n')
      prmfile.write(' save2Dbin    = .true.\n\n')

      prmfile.write(' nsave3D      = 20,\n')
      prmfile.write(' nsave3Dstart	= 20,\n')
      prmfile.write(' nsave3Dend   = 999999999,\n')
      prmfile.write(' save3Dbin    = .true.\n\n')

      prmfile.write(' nstatmom      = 20\n')
      prmfile.write(' nstatmomstart = 9999999\n')
      prmfile.write(' nstatmomend   = 99960480\n')
      prmfile.write(' savemombin    = .true.\n\n')

      prmfile.write(' perturb_type = 5, \n\n')

      prmfile.write(' dotracerperiodic = .true.,\n\n')

      prmfile.write(' massemis = '+str(massemis)+',\n')
      prmfile.write(' vg = '+str(vg)+', \n')
      prmfile.write(' plume_y = 1000.,\n')     

      prmfile.write(' /\n')         

      prmfile.close()


      # --- write the grd file ---
      print 'writing grd file...'
      outputfile =  open(rundir+'/RUNS/'+casename +'/grd','w')
      outputfile.write('20\n60\n100')
      outputfile.close()

      # --- copy in lst file ---
      # Maybe change this to a series of write statements later.
      # It's 188 lines, so I don't really want that all in here.
      print 'copying in lst file...'
      #os.system('cp lst '+rundir+'/RUNS/'+casename+'/')
      os.system('cp '+rundir+'/RUNDATA/lst '+rundir+'/RUNS/'+casename+'/') 

      # --------------------------------------------------------------------
      # 3) Queue the model to be run. This problem was solved for
      # StrawberryFields, (this was in my queueSAM.py script. Now I must
      # re-solve it for running on the ACE-Net machines, using their
      # queueing system.
      # --------------------------------------------------------------------

               
   caseid+=1
