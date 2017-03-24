#!/usr/bin/env python
#===============================================================================

#===============================================================================

import sys, os
import string
import re
import glob
import matplotlib.pylab as P
import numpy as N
import time
import netCDF4 as ncdf
from pyproj import Proj
from optparse import OptionParser
from tables import *
from netcdftime import utime
from datetime import datetime as py_datetime
from scipy import interpolate
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid import AxesGrid
from mpl_toolkits.axes_grid.inset_locator import inset_axes

_debug               = True
_verbose             = True
_missing             = -999.
version_string       = "pyDART_file_version_3.0"
checked_file_version = False

#=========================================================================================
# PARAMETERS FOR MRMS

dbz_zeros         = True            # True: missing < dBZ < dbz_min=>0, False: set to missing (creates halo)
dbz_missing_zeros = True            # True => fill _missing with zeros, False => OFF
dbz_min           = 20.             # Reflectivity values below this are _missing or set to dbz_clear_air
dbz_max           = 70.             # Set reflectivity values above this value to dbz_max
dbz_thin          = 4               # thin data by skipping this in x&y (1 ==> no thinning, usually 3 or 6)
dbz_thin_z        = 2               # Vertically thin data by this factor (usually 2) 
dbz_clear_thin    = 2               # Thin zeros additionally by skipping this in x&y (1 ==> no thinning)
dbz_clear_air     = 0.0             # value for clear air reflectivity (usually 0 or -10)
dbz_clear_hgts    = (True, 3000., 7000.)
dbz_stdev         = 7.5
lat_bound_box     = (33.0,37.0)     # max and min latitudes to output
lon_bound_box     = (-100.,-96.0)   # max and min longitudes to output
hgt_bound_box     = (0.0,10000.0)   # max and min heights to output

#==========================================================================================
# PARAMETERS FOR MAP PROJECTIONS

map_projection     = 'lcc'  # valid projects are lambert conformal ['lcc'], and latlon ['latlon']
truelat1, truelat2 = 30.0, 60.0

#==========================================================================================
# Some common radar locations that can be handy

radar_LookUp = { 
                 "KOUN":          (35.2358, -97.4622, 370.),
                 "MPAR":          (35.2358, -97.4622, 370.),
                 "KTLX":          (35.333,  -97.278,  415.),
                 "KFDR":          (34.3520, -98.9839, 383.),
                 "KVNX":          (36.7408, -98.1275, 369.),
                 "SR2_30May2004": (35.4976, -98.4654, -999.)
               }

# default radar location

_radar_loc =  radar_LookUp["MPAR"][:]

#==========================================================================================
# Dont touch below this file

day_utime   = utime("days since 1601-01-01 00:00:00")
sec_utime   = utime("seconds since 1970-01-01 00:00:00")
time_format = "%Y-%m-%d_%H:%M:%S"

default_start = (1984, 10, 12, 7, 35, 30)
default_end   = (2020, 10, 12, 7, 35, 30)

hscale        = 1./1000.  # horizontal scale for plotting

#===============================================================================
def nice_mxmnintvl( dmin, dmax, outside=False, max_steps=20, cint=None, sym=False):

    """ Description: Given min and max values of a data domain and the maximum
                     number of steps desired, determines "nice" values of
                     for endpoints and spacing to create a series of steps
                     through the data domainp. A flag controls whether the max
                     and min are inside or outside the data range.

        In Args: float   dmin             the minimum value of the domain
                 float   dmax       the maximum value of the domain
                 int     max_steps  the maximum number of steps desired
                 logical outside    controls whether return min/max fall just
                                    outside or just inside the data domainp.
                     if outside:
                         min_out <= min < min_out + step_size
                                         max_out >= max > max_out - step_size
                     if inside:
                         min_out >= min > min_out - step_size
                                         max_out <= max < max_out + step_size

                 float    cint      if specified, the contour interval is set
                                    to this, and the max/min bounds, based on
                                    "outside" are returned.

                 logical  sym       if True, set the max/min bounds to be anti-symmetric.


        Out Args: min_out     a "nice" minimum value
                  max_out     a "nice" maximum value
                  step_size   a step value such that
                                     (where n is an integer < max_steps):
                                      min_out + n * step_size == max_out
                                      with no remainder

        If max==min, or a contour interval cannot be computed, returns "None"

        Algorithm mimics the NCAR NCL lib "nice_mxmnintvl"; code adapted from
        "nicevals.c" however, added the optional "cint" arg to facilitate user
        specified specific interval.

        Lou Wicker, August 2009 """

    table = N.array([1.0,2.0,2.5,4.0,5.0,10.0,20.0,25.0,40.0,50.0,100.0,200.0,
                      250.0,400.0,500.0])

    if nearlyequal(dmax,dmin):
        return None, None, N.zero(1)

    # Help people like me who can never remember - flip max/min if inputted reversed
    if dmax < dmin:
        amax = dmin
        amin = dmax
    else:
        amax = dmax
        amin = dmin

    if sym:
        smax = max(amax.max(), amin.min())
        amax = smax
        amin = -smax
        
    d = 10.0**(N.floor(N.log10(amax - amin)) - 2.0)
    if cint == None or cint == 0.0:
        t = table * d
    else:
        t = cint
    if outside:
        am1 = N.floor(amin/t) * t
        ax1 = N.ceil(amax/t)  * t
        cints = (ax1 - am1) / t
    else:
        am1 = N.ceil(amin/t) * t
        ax1 = N.floor(amax/t)  * t
        cints = (ax1 - am1) / t

    # DEBUG LINE BELOW
    #print am1
    #print ax1
    #print cints

    if cint == None or cint == 0.0:
        try:
            index = N.where(cints < max_steps)[0][0]
            return am1[index], ax1[index], N.linspace(am1[index], ax1[index], cints[index])
        except IndexError:
            return None, None, N.zero(1)
    else:
        return am1, ax1, N.arange(am1, ax1+cint, cint)
#===============================================================================        
def nearlyequal(a, b, sig_digit=None):

    """ Measures the equality (for two floats), in unit of decimal significant
        figures.  If no sigificant digit is specified, default is 7 digits. """

    if sig_digit == None or sig_digit > 7:
        sig_digit = 7
    if a == b:
        return True
    difference = abs(a - b)
    avg = (a + b)/2

    return N.log10(avg / difference) >= sig_digit

#===============================================================================        
def beam_elv(sfc_range, z):
########################################################################
#
#     PURPOSE:
#
#     Calculate the elevation angle (elvang) and the along
#     ray-path distance (range) of a radar beam
#     crossing through the given height and along-ground
#     distance.
#
#     This method assumes dn/dh is constant such that the
#     beam curves with a radius of 4/3 of the earth's radius.
#     This is dervied from Eq. 2.28 of Doviak and Zrnic',
#     Doppler Radar and Weather Observations, 1st Ed.
#
########################################################################
#
#     AUTHOR: Keith Brewster
#     10/10/95
#
#     MODIFICATION HISTORY: adapted to python by Lou Wicker (thanks Keith)
#
########################################################################
#
#     INPUT:
#       sfc_range:    Distance (meters) along ground from radar 
#       z        :    Height above radar
#
#     OUTPUT
#       elvang   Elevation angle (degrees) of radar beam
#
########################################################################
    eradius=6371000.
    frthrde=(4.*eradius/3.)
    eighthre=(8.*eradius/3.)
    fthsq=(frthrde*frthrde)

    if sfc_range > 0.0:
        hgtdb = frthrde + z
        rngdb = sfc_range/frthrde

        elvrad = N.arctan((hgtdb*N.cos(rngdb) - frthrde)/(hgtdb * N.sin(rngdb)))

        return N.rad2deg(elvrad)

    else:

        return -999.

#===============================================================================
def dxy_2_dll(x, y, lat1, lon1, degrees=True, proj = map_projection):

  """dxy_2_dll returns the approximate lat/lon between an x,y coordinate and
     a reference lat/lon point.  

     Valid projections: Lambert Conformal 
                        Lat - Lon

     INPUTS:  x,y in meters, lat1, lon1 in radians, or if degrees == True, 
              then degrees (default value)

     if x > 0, lon > lon1

     if y > 0, lat > lat1

     OUTPUTS:  lat, lon in radians, or if degrees == True, degrees
               i.e., the input and output units for lat/lon are held the same.
  """

  if degrees:
    rlon1 = N.deg2rad(lon1)
    rlat1 = N.deg2rad(lat1)
  else:
    rlon1 = lon1
    rlat1 = lat1

# Simple Lat-Lon grid

  if proj == 'latlon':
    rearth = 1000.0 * 6367.0
    rlat2  = rlat1 + y / rearth
    lon    = N.rad2deg(rlon1 + x / ( rearth * N.cos(0.5*(rlat1+rlat1)) ) )
    lat    = N.rad2deg(rlat2)

# Lambert Conformal

  if proj == 'lcc':
    p1 = Proj(proj='lcc', ellps='WGS84', datum='WGS84', lat_1=truelat1, lat_2=truelat2, lat_0=lat1, lon_0=lon1)
    lon, lat = p1(x, y, inverse = True)

  if degrees == False:
    return N.deg2rad(lat), N.deg2rad(lon)
  else:
    return lat, lon

#===============================================================================
def dll_2_dxy(lat1, lat2, lon1, lon2, degrees=False, azimuth=False, proj = map_projection):

  """dll_2_dxy returns the approximate distance in meters between two lat/lon pairs

     Valid projections: Lambert Conformal 
                        Lat - Lon

     INPUTS: Two (lat,lon) pairs in radians, or if degrees==True, degrees (default)

     if lon2 > lon1:  x > 0

     if lat2 > lat1:  y > 0

     OUTPUTS:  DX, DY in meters

     Azimuth formula from http://www.movable-type.co.uk/scripts/latlong.html
  """

  if degrees:
    rlon1 = N.deg2rad(lon1)
    rlon2 = N.deg2rad(lon2)
    rlat1 = N.deg2rad(lat1)
    rlat2 = N.deg2rad(lat2)
  else:
    rlon1 = lon1
    rlon2 = lon2
    rlat1 = lat1
    rlat2 = lat2

# Simple Lat-Lon grid

  if proj == 'latlon':
    rearth  = 1000.0 * 6367.0
    x       = rearth * N.cos(0.5*(rlat1+rlat2)) * (rlon2-rlon1)
    y       = rearth * (rlat2-rlat1)

# Lambert Conformal

  if proj == 'lcc':
    p1 = Proj(proj='lcc', ellps='WGS84', datum='WGS84', lat_1=truelat1, lat_2=truelat2, lat_0=lat1, lon_0=lon1)
    x, y = p1(lon2, lat2, errchk = True)

  if azimuth:
    ay = N.sin(rlon2-rlon1)*N.cos(rlat2)
    ax = N.cos(rlat1)*N.sin(rlat2)-N.sin(rlat1)*N.cos(rlat2)*N.cos(rlon2-rlon1)
    az = N.degrees(N.arctan2(ay,ax))
    return x, y, az

  return x, y

#===================================================================================================
def cressman(x, y, obs, x0, y0, roi, missing=_missing):
  """ Returns a data value for the point
      Arguments: x/y/obs:  1D arrays of location
                 x0, y0:   point to analyze to
                 roi:      radius of influence
                 missing:  value to assign if no data, default = 0.0
      This routine can also be implemented in FORTRAN much quicker...
  """
# Create distance array

  dis = N.sqrt( (x-x0)**2 + (y-y0)**2 )

# Cut the size o n by choosing a smart threshold (2 km)

  indices = N.where(dis <= roi)

# Check to see if there are any data ponts that are within search radius

  size = N.size(indices)

  if size != 0:
    R2    = roi**2.0
    w_sum = 0.0
    top   = 0.0
    print "CRESSMAN param:  ", size, x, y, obs, x0, y0

# go thru values w/in radius: calculate weight, mult by value and sum also sum weights

    for n, value in enumerate(dis[indices]):
      rk2 = value**2.0
      wk = (R2-rk2) / (R2+rk2)
      top = top + wk*obs[n]
      w_sum = w_sum + wk
#       print n, value, data[n], R2, w_sum, top

    if w_sum < 0.01:
      return missing
    else:
      return top/w_sum

  else:
    return missing

#===============================================================================
def chk_pyDart_version(h5file, verbose = True):
    """Checks the H5 pyDart file to see whether the file format is compatible"""
    
    global checked_file_version
    
    if checked_file_version and h5file.title != version_string:
        print
        print "\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/"
        print " !!! pyDART file is wrong version !!!"
        print "pyDart sofware expects version:  ", version_string
        print "pyDart file contains version:    ", h5file.title
        print "pyDart software is exiting"
        print "/\/\/\/\/\/\/\/\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/"
        print
        sys.exit(0)
    else:
        if not checked_file_version and verbose == True:
            print
            print "pyDART file is version:  ",h5file.title
            print "/\/\/\/\/\/\/\/\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/"
            print
            return

#===============================================================================
def open_pyDart_file(filename, return_root=False, verbose = None, append=False):
    
    """Open_pyDart_file opens the pyTables pyDart file, checks a few things and
       returns a the file handle and pyDart table object for use.
    """

# Open DART PyTables file
    
    if append:
        h5file = open_file(filename, mode = "a")
    else:
        h5file = open_file(filename, mode = "r")
    
    chk_pyDart_version(h5file,verbose)
    
    root   = h5file.root

# observations are a table in obs.observations
    
    if return_root:
        return h5file, root, table
    else:
        return h5file, root.obs.observations


#===================================================================================================
def ObType_LookUp(name,DART_name=False,Print_Table=False):
      """ObType_LookUp returns the DART kind number for an input variable type.  There seems
         to be several ways observations names are written in the observation inputs
         and output files, e.g., in the DART ascii files and in the ***.obs.nc files,
         so this function is designed to handle the variety of cases and return the
         integer corresponding to the DART definition.
      
         Exampled:   REFLECTIVITY is sometimes stored as REFL
                     T_2M         is sometimes stored as TEMP2m
                     TD_2M        is sometimes stored as DEWPT2m
      
         If you come across a variable that is not defined, you can add it to the lookup
         table (dictionary) below - just make sure you know the official DART definition
         
         You can add any unique variable name or reference on the left column, and then use
         the pyDART INTERNAL definition on the right - so that you can refer to more than
         one type of data in different ways internally in this code
      
         Usage:  variable_kind = ObType_Lookup(variable_name)   where type(variable_name)=str
      
         If you need the return the actual DART_name as well, set the input flag to be True"""

# Create local dictionary for observation kind definition - these can include user abbreviations

#                      user's observation type            kind   DART official name

      Look_Up_Table={ "DOPPLER_VELOCITY":                 [11,   "DOPPLER_RADIAL_VELOCITY"] ,
                      "DOPPLER_RADIAL_VELOCITY":          [11,   "DOPPLER_RADIAL_VELOCITY"] ,
                      "REFLECTIVITY":                     [12,   "RADAR_REFLECTIVITY"],
                      "RADAR_REFLECTIVITY":               [12,   "RADAR_REFLECTIVITY"],
                      "RADAR_CLEARAIR_REFLECTIVITY":      [13,   "RADAR_CLEARAIR_REFLECTIVITY"],  
                      "CLEARAIR_REFLECTIVITY":            [13,   "RADAR_CLEARAIR_REFLECTIVITY"],
                      "DIFFERENTIAL_REFLECTIVITY":        [300,  "DIFFERENTIAL_REFLECTIVITY"],
                      "SPECIFIC_DIFFERENTIAL_PHASE":      [301,  "SPECIFIC_DIFFERENTIAL_PHASE"],
                      "METAR_U_10_METER_WIND":            [1,    "METAR_U_10_METER_WIND"],
                      "METAR_V_10_METER_WIND":            [2,    "METAR_V_10_METER_WIND"],
                      "METAR_TEMPERATURE_2_METER":        [4,    "METAR_TEMPERATURE_2_METER"],
                      "METAR_DEWPOINT_2_METER":           [9,    "METAR_DEWPOINT_2_METER"],
                      "METAR_SPECIFIC_HUMIDITY_2_METER":  [5,    "METAR_SPECIFIC_HUMIDITY_2_METER"],
                      "VR":                               [11,   "DOPPLER_RADIAL_VELOCITY"],
                      "DBZ":                              [12,   "RADAR_REFLECTIVITY"],
                      "0DBZ":                             [13,   "RADAR_CLEARAIR_REFLECTIVITY"],
                      "ZDR":                              [300,  "DIFFERENTIAL_REFLECTIVITY"],
                      "KDP":                              [301,  "SPECIFIC_DIFFERENTIAL_PHASE"],
                      "U10M":                             [1,    "METAR_U_10_METER_WIND"],
                      "V10M":                             [2,    "METAR_V_10_METER_WIND"],
                      "T2M":                              [4,    "METAR_TEMPERATURE_2_METER"],
                      "TD2M":                             [9,    "METAR_DEWPOINT_2_METER"],
                      "H2M":                              [5,    "METAR_SPECIFIC_HUMIDITY_2_METER"],
                      "U_10M":                            [1,    "METAR_U_10_METER_WIND"],
                      "V_10M":                            [2,    "METAR_V_10_METER_WIND"],
                      "T_2M":                             [4,    "METAR_TEMPERATURE_2_METER"],
                      "TD_2M":                            [9,    "DEW_POINT_2_METER"],
                      "H_2M":                             [5,    "METAR_SPECIFIC_HUMIDITY_2_METER"],
                      "TEMP2M":                           [4,    "METAR_TEMPERATURE_2_METER"],
                      "DEWPT2M":                          [9,    "DEW_POINT_2_METER"],
                      "REFL":                             [12,   "RADAR_REFLECTIVITY"],
                      "FLASH_RATE_2D":                    [2014, "FLASH_RATE_2D"],
                      "METAR_ALTIMETER":                  [71,   "METAR_ALTIMETER"],
                      "LAND_SFC_ALTIMETER":               [70,   "LAND_SFC_ALTIMETER"],
                      "LAND_SFC_DEWPOINT":                [58,   "LAND_SFC_DEWPOINT"],
                      "LAND_SFC_TEMPERATURE":             [25,   "LAND_SFC_TEMPERATURE"],
                      "LAND_SFC_U_WIND_COMPONENT":        [23,   "LAND_SFC_U_WIND_COMPONENT"],
                      "LAND_SFC_V_WIND_COMPONENT":        [24,   "LAND_SFC_V_WIND_COMPONENT"],
                      "GOES_CWP_PATH":                    [80,   "GOES_CWP_PATH"],
                      "GOES_LWP_PATH":                    [119,   "GOES_LWP_PATH"],
                      "GOES_IWP_PATH":                    [120,   "GOES_IWP_PATH"],
                      "GOES_CWP_ZERO":                    [121,   "GOES_IWP_ZERO"]
                    }
      
      if Print_Table:
            print
            print "VALID INPUT VARIABLE NAME              KIND  DART NAME"
            print "=========================================================================="
            for key in Look_Up_Table.keys():
                  print "%35s    %3d    %s" % (key, Look_Up_Table[key][0], Look_Up_Table[key][1])
            return
      
      name2 = name.upper().strip()
      if Look_Up_Table.has_key(name2):
            if DART_name == True:
                  return Look_Up_Table[name2][0], Look_Up_Table[name2][1]
            else:
                 return Look_Up_Table[name2][0]
      else:
            print "ObType_LookUp cannot find variable:  ", name, name2
            raise SystemExit

#===============================================================================
def mergeTables(table_new, tables, addindex=True):

    print "mergeTable called:  New table:  ",table_new

    if len(tables) == 2 and tables[0] == "D":
        tables = glob.glob("%s/*.h5" % tables[1])

    print "mergeTable called:  Reading from:  ",tables, len(list(tables))

    if len(tables) == 1 or type(tables) == type('str'):
        print "Only one table for merging supplied, simply doing a copy..."
        h5file1, table1 = open_pyDart_file(tables)
        h5file1.copy_file(table_new, overwrite=True)
        h5file1.close()
        print "Finished copying %s into %s" % (tables, table_new)
    else:
        h5file1, table1 = open_pyDart_file(tables[0])
        print "Creating new table to copy into...."
        h5file1.copy_file(table_new, overwrite=True)
        print("Finished copying %s into %s, has a length: %i" % (tables[0], table1, table1.nrows))

        h5file1.close()   
        h5file1, table1 = open_pyDart_file(table_new, append=True)

        for file in tables[1:]: 
            print "Now copying from table:  ", file
            h5file2, table2 = open_pyDart_file(file)

            rows = table2.read()
            table1.append(rows)
            table1.flush()

            h5file2.close()
            print("Finished copying %s into %s\n New table has a length: %i\n" % (tables[0], table1, table1.nrows))

        print "Finished appending all table rows...."
        print "New table:    ", table1
        
        indexrows = table1.cols.utime.create_csindex()
        
        group_header = h5file1.root.header
        group_header.attributes.cols.last[0]        = table1.nrows
        group_header.attributes.cols.max_num_obs[0] = table1.nrows
        group_header.attributes.cols.num_obs[0]     = table1.nrows
        
        h5file1.close()

        return

#===============================================================================
#
def sortTable(filename, overwrite=True):

    if overwrite:
       cmd = ("cp %s %s_unsorted.h5" % (filename, filename[:-3]))
       print("\n sortTable is running command:  %s" % cmd)
       os.system(cmd)
       
    cmd = ("ptrepack -o --overwrite-nodes --keep-source-filters %s:/ sorted.h5:/" % filename)
    print("\n sortTable is running command:  %s" % cmd)
    os.system(cmd)
    
    cmd = ("ptrepack --sortby utime --overwrite-nodes --keep-source-filters %s:/obs/observations sorted.h5:/obs/observations" % \
           filename)
    print("\n sortTable is running command:  %s" % cmd)
    os.system(cmd)
    
    if overwrite:
         cmd = ("mv sorted.h5 %s" % filename)
         print("\n sortTable is running command:  %s" % cmd)
         os.system(cmd)
                   
    return

#===============================================================================
class pyDART():

    def __init__(self, verbose=True, debug=True):
        self.hdf5    = None
        self.ascii   = None
        self.index   = None
        self.verbose = verbose
        self.debug   = debug

#-------------------------------------------------------------------------------
# File:  creates a standard set of filenames for pyDART ascii and hdf5 files
#-------------------------------------------------------------------------------
    
    def file(self,filename=None):
        
        if filename != None:
            
            if filename[-4:] == ".out":
                self.hdf5  = filename[:-4] + ".h5"
                self.ascii = filename
                return
            
            if filename[-3:] == ".h5":
                self.hdf5  = filename
                self.ascii = filename[:-3] + ".out"
                return
            
            if filename[-5:] == ".hdf5":
                self.hdf5  = filename
                self.ascii = filename[:-5] + ".out"
                return
        
        else:
            print "pyDART.init:  No file is supplied, please add one to self.hdf5 or self.ascii"
        
        return
#-------------------------------------------------------------------------------
# Search:
#-------------------------------------------------------------------------------
    
    def search(self, variable=None, start=None, end=None, condition=None, loc=None, selfdata=False, tablereturn=None):

# Construct a variable to search table

        if variable == None:
            print "pyDART.search:  no variable supplied....all vars dumped"
        if start == None and end == None:
            print "pyDART.search:  no start or end search time supplied.... using default"
            start = default_start
            end = default_end
        
        if start != None and type(start) == tuple:
            self.start = py_datetime(start[0],start[1],start[2],start[3],start[4],start[5])
            if self.verbose:
                if self.debug:  print "pyDART.search:  Start datetime:  ",self.start
        
        if end != None and type(end) == tuple:
            self.end = py_datetime(end[0],end[1],end[2],end[3],end[4],end[5])
            if self.verbose:
                if self.debug:  print "pyDART.search:  End datetime:  ", self.end

# Build a list of search conditions

        cond = []
        
        if variable != None:
            cond.append( "( kind == " + str(ObType_LookUp(variable)) + " )" )
        
        if start != None:
            utime_start = sec_utime.date2num(self.start)
            cond.append( "(" + str(utime_start) + " <= utime)" )
        
        if end != None:
            utime_end   = sec_utime.date2num(self.end)
            cond.append( "(utime < " + str(utime_end) + ")" )
        
        if condition != None:
            cond.append(condition)
        
        if loc != None:
            for item in loc:
                cond.append(item)

# Open DART PyTables file
        
        h5file, table = open_pyDart_file(self.hdf5, verbose=self.verbose)

# Determine when you want observations from
        
        if self.verbose:
            print "PyDART SEARCH START TIME  %s  /  UTIME_START: %s" % (self.start, utime_start)
            print "PyDART SEARCH END   TIME  %s  /  UTIME_END:   %s" % (self.end, utime_end)
            print "PyDART converted utimes: ", sec_utime.num2date(utime_start)

# Create string for search
        
        if len(cond) != 0:
            
            search_string = cond[0]
            for x in cond[1:]:
                search_string = search_string + " & " + x
            
            if self.verbose:
                print
                print "PyDART SEARCH CONDITION IS:  ", search_string
                print
            
            self.index = table.get_where_list(search_string)    # Do the search

            if len(self.index) == 0:  self.index = None
            
            if tablereturn != None:

                # make a blank table to put results of search in
                filter_spec = Filters(complevel=5, complib="zlib", shuffle=1, fletcher32=0)
                h5file_sub = open_file(tablereturn, mode = "w", title = version_string, filters=filter_spec)
                group_ob_kinds = h5file_sub.create_group("/", 'obs', 'Obs for DART file')
        
                table_ob_kinds = h5file_sub.create_table(group_ob_kinds, 'kinds', DART_ob_kinds, 'Observation Descriptions')
                group_header = h5file_sub.create_group("/", 'header', 'Header Information for DART file')
                table_header = h5file_sub.create_table(group_header, 'attributes', DART_header, 'Attributes of the observational file')
        
                root         = h5file_sub.root
                group_obs    = root.obs
                group_header = root.header
        
                table_obs = h5file_sub.create_table(group_obs, 'observations', DART_obs, 'Observations from DART file')
                
                # do search and put results in table
                table.append_where(table_obs, search_string)
                
                table_obs.flush()
                
                
                h5file_sub.close()
        
        else:
            if self.verbose:
                print "pyDART.SEARCH:  NO SEARCH CONDITION CREATED....", cond
        
        h5file.close()
        
        return

#-------------------------------------------------------------------------------
# A quick routine to grid pyDART data
    
#-------------------------------------------------------------------------------
# A quick routine to grid pyDART data
    
    def grid(self, xo, yo, zo, dx=None, dy=None):

# We assume that the data has been thinned to remove duplicates, or close enough....
        
# Always hard to create a reasonable guess at a grid.  Find middle of data, and work back
#        from there.  Using geometric, not arithmatic, mean for domain size.  Use 80% of it.
        
        xmean = xo.mean()
        ymean = yo.mean()
        del_y = yo.max() - yo.min()
        del_x = xo.max() - xo.min()
        xmin  = xmean - 1.2*(del_y*del_x / (del_y+del_x))
        ymin  = ymean - 1.2*(del_y*del_x / (del_y+del_x))
        xmax  = xmean + 1.2*(del_y*del_x / (del_y+del_x))
        ymax  = ymean + 1.2*(del_y*del_x / (del_y+del_x))

        if dx == None:
            dx = del_x / 50.
        if dy == None:
            dy = del_y / 50.

        xi = N.arange(xmin, xmax, dx)
        yi = N.arange(ymin, ymax, dy)

# Radial basis function method - not so good.
#       rbf_z = Rbf(xo, yo, zo, function='thin_plate', epsilon=2)
#       xx, yy = N.meshgrid(xi, yi)
#       zi     = rbf_z(xx,yy)
#
#       zi = interpolate.griddata((xo, yo), zo, (xi[None,:], yi[:,None]), method='cubic')
        zi = P.mlab.griddata(xo,yo,zo,xi,yi,interp='linear')
        
        if self.verbose:
            print
            print 'X(nx) min/max of grid region:     ',xi.size, xi.min(), xi.max()
            print 'Y(ny) min/max of grid region:     ',yi.size, yi.min(), yi.max()
            print 'Z(nx,ny) min/max of grid region:  ',zi.shape, zi.min(), zi.max()
        
        return xi, yi, zi

#-------------------------------------------------------------------------------
# A scatter plot routine
    
    def scatter(self, variable=None, convertLatLon=False):

# Return some data specified (hopefully) by user search
# Later, we might add an error check here to make sure the data are all the same kind
        
        data = self.get_data()
        time = data['date'][0]

# Need to remove all locations that are identical - do that by creating a complex number to
#      comprised of (x + j * y) which can then be sorted for multiple identical entries.
#      We then remove those entries and grid them using Jeff Whitaker's nat grid which has been
#      added to the ENTHOUGHT/numpy libs
#      (see http://sourceforge.net/project/showfiles.php?group_id=80706&package_id=142792)
        
        xy = data['lon'].copy() + (0.+1.j)*data['lat'].copy()
        xy_unique, index = N.unique(xy, return_index=True)
        
        print
        print "Number of points found from initial search:      %d" % (xy.size)
        print "Number of points which are unique (no overlap):  %d" % (xy_unique.size)

        x   = N.real(xy_unique)
        y   = N.imag(xy_unique)
        z   = data['value'][index[:]]

        if convertLatLon:
           x = N.rad2deg(x)
           y = N.rad2deg(y)
    
        x   =  N.where(x > 180., x - 360., x)  # make west longitudes...
        z   =  N.array(z, dtype="float32")

        el = N.round(data['elevation'].mean(),2)

        if self.verbose:
            print
            print 'X(nx) min/max of grid region:     ',x.size, x.min(), x.max()
            print 'Y(ny) min/max of grid region:     ',y.size, y.min(), y.max()
            print 'Z(nx,ny) min/max of grid region:  ',z.shape, z.min(), z.max()


        xmax = x.max()
        ymax = y.max()
        xmin = x.min()
        ymin = y.min()
    
# color map
        Zmap  = [  (0.000,1.000,1.000), \
                    (0.118,0.566,1.000), \
                    (0.000,0.000,0.804), \
                    (0.482,0.988,0.000), \
                    (0.000,0.933,0.000), \
                    (0.000,0.545,0.000), \
                    (1.000,1.000,0.000), \
                    (0.726,0.525,0.043), \
                    (1.000,0.549,0.000), \
                    (1.000,0.000,0.000), \
                    (0.804,0.000,0.000), \
                    (0.549,0.000,0.000), \
                    (0.933,0.070,0.537), \
                    (0.604,0.196,0.078) ]
    
        fig = P.figure(figsize = (16,10))
        ax = fig.add_subplot(111, aspect='equal')
        
        sw_lon = x.min()
        ne_lon = x.max()
        sw_lat = y.min()
        ne_lat = y.max()
        
        map = Basemap(projection='lcc', llcrnrlon=sw_lon,llcrnrlat=sw_lat,urcrnrlon=ne_lon,urcrnrlat=ne_lat, \
              lat_0=0.5*(ne_lat+sw_lat), lon_0=0.5*(ne_lon+sw_lon), resolution='l',area_thresh=1.,suppress_ticks=True)
        
        map.drawstates(linewidth=2.0, color='k')
        map.drawcounties(linewidth=0.5, linestyle='solid', color='k')
        
        xx, yy = map(x, y)
        
        if variable.upper() == "DBZ" or variable.upper() == "REFLECTIVITY":
            plt = map.scatter(xx, yy, c=z, vmin=0, vmax=70, edgecolors='none')
            cbar = map.colorbar(plt,location='bottom',pad=0.40)
            cbar.ax.set_title("dBZ")
            P.title("Reflectivity at   %s  EL:  %5.2f" % (time, el))
            
        elif variable.upper() == "0DBZ" or variable.upper() == "CLEAR_AIR_REFLECTIVITY":
            plt = map.scatter(xx,yy,c=z, vmin=0, vmax=1, edgecolors='none')
            cbar = map.colorbar(plt,location='bottom',pad=0.40)
            cbar.ax.set_title("dBZ")
            P.title("0-Reflectivity at   %s  EL:  %f" % (time, el))
            
        else:
            amin1, amax1, clevels = nice_mxmnintvl(z.min(), z.max())
            plt  = map.scatter(xx,yy,c=z, vmin=amin1, vmax=amax1, edgecolors='none',cmap='jet')
            cbar = map.colorbar(plt,location='bottom',pad=0.40)
            cbar.ax.set_title(variable)
            P.title("%s at   %s  EL:  %f" % (variable, time, el))

        P.savefig('DART_scatter.png')
        P.show()
    
        return

# end modifications RLT 20110113
#-------------------------------------------------------------------------------
# A quick and dirty plotting routine to make sure the pyDART is reasonable

    def plot(self, variable=None, convertLatLon=False, savefig=None):

# define a quick and dirty colormap for reflectivity
        
        cmap  = [  (0.000,1.000,1.000), \
                (0.118,0.566,1.000), \
                (0.000,0.000,0.804), \
                (0.482,0.988,0.000), \
                (0.000,0.933,0.000), \
                (0.000,0.545,0.000), \
                (1.000,1.000,0.000), \
                (0.726,0.525,0.043), \
                (1.000,0.549,0.000), \
                (1.000,0.000,0.000), \
                (0.804,0.000,0.000), \
                (0.549,0.000,0.000), \
                (0.933,0.070,0.537), \
                (0.604,0.196,0.078) ]

# Return some data specified (hopefully) by user search
# Later, we might add an error check here to make sure the data are all the same kind
        
        data = self.get_data()
        time = data['date'][0]

# Need to remove all locations that are identical - do that by creating a complex number to
#      comprised of (x + j * y) which can then be sorted for multiple identical entries.
        
        xy = data['lon'].copy() + (0.+1.j)*data['lat'].copy()
        xy_unique, index = N.unique(xy, return_index=True)
        
        print "\nNumber of points found from initial search:      %d" % (xy.size)
        print "Number of points which are unique (no overlap):  %d\n" % (xy_unique.size)

        x = N.real(xy_unique)
        y = N.imag(xy_unique)
        z = data['value'][index[:]]

        if convertLatLon:
           x = N.rad2deg(x)
           y = N.rad2deg(y)
    
        x =  N.where(x > 180., x - 360., x)  # make west longitudes...
        z =  N.array(z, dtype="float32")
        
        el = N.round(data['elevation'].mean(),2)

        print "Creating grid for elevation:  %f" % data['elevation'].mean()
        
# Now, call grid to grid some data
        
        xi, yi, zi = self.grid(x, y, z, dx=0.02, dy=0.02)

        fig = P.figure(figsize = (12,10))

        ax = fig.add_subplot(111)
        
        sw_lon = xi.min()
        ne_lon = xi.max()
        sw_lat = yi.min()
        ne_lat = yi.max()
        
        map = Basemap(projection='lcc', llcrnrlon=sw_lon,llcrnrlat=sw_lat,urcrnrlon=ne_lon,urcrnrlat=ne_lat, \
              lat_0=0.5*(ne_lat+sw_lat), lon_0=0.5*(ne_lon+sw_lon), resolution='i',area_thresh=1.)
        
        lon2d, lat2d = N.meshgrid(xi, yi)        
        xx, yy = map(lon2d, lat2d)
        
        map.drawcounties(linewidth=0.5, linestyle='solid', color='gray')
        map.drawstates(linewidth=2.0, color='k')

        if variable.upper() == "DBZ" or variable.upper() == "REFLECTIVITY":
            clevels = N.arange(0,75,5)
            zi   = N.ma.masked_less( zi, 10.0 )
            plt  = map.contourf(xx, yy, zi, clevels, colors=cmap)
            cbar = map.colorbar(plt,location='bottom',pad=0.40)
            cbar.ax.set_title("dBZ")
            P.title("Reflectivity at   %s  EL:  %5.2f" % (time, el))

# Plot 0's values as points.

            r_mask = (zi.mask == True)
            map.scatter(xx[r_mask], yy[r_mask], c='k', s = 1., alpha=0.5)

        elif variable.upper() == "0DBZ" or variable.upper() == "CLEAR_AIR_REFLECTIVITY":
            clevels = N.arange(0,10,5)
            plt = map.contourf(xx,yy,zi, clevels, colors=cmap)
            cbar = map.colorbar(plt,location='bottom',pad=0.40)
            cbar.ax.set_title("Clear Air Reflecivity")
            P.title("0-Reflectivity at   %s  EL:  %5.2f" % (time, el))
            
        else:
            amin1, amax1, clevels = nice_mxmnintvl(zi.min(), zi.max())
            plt = map.contourf(xx, yy, zi, clevels)
            cbar = map.colorbar(plt,location='bottom',pad=0.40)
            cbar.ax.set_title(variable)
            P.title("%s at   %s  EL:  %5.2f" % (variable, time, el))

        P.savefig("DART_contour.png")

        P.show()
        
        return

#-------------------------------------------------------------------------------
    
    def get_data(self, variable=None, all=False):
        
        if self.index == None and all ==False:
            if self.verbose: print "pyDART.get_data:  No search indices supplied, returning all rows"

# Open DART PyTables file
        
        h5file, table = open_pyDart_file(self.hdf5, verbose=self.verbose)
        
        if self.index == None and all == True:
            self.index = N.arange(table.nrows)
            if self.verbose: print "pyDART.get_data, return all rows of table!"
        
        if self.debug:  start = time.clock()
        
        data = table.read_coordinates(self.index)

        if self.debug:  print "pyDart.get_data  Execution time for read_coordinates method:  ", time.clock() - start
         
        h5file.close()
        
        return data

#-------------------------------------------------------------------------------
    
    def list(self,variable=None,dumplength=None):
        
        if variable == None:
            if self.verbose: print "pyDART.list:  No variable supplied, listing all variables in ", self.hdf5,"/obs/observations"
            variable = "all"

# Open DART PyTables file
        
        h5file, table = open_pyDart_file(self.hdf5, verbose=self.verbose)
        
        if type(variable) == type("str"):       # Determine the user's request data type
            if variable == "all":               # Retrieve all the columns from the row
                for name in table.colnames:
                    print name
                h5file.close()
                return
            else:                                     # User requests specific variables, search for DART variable index
                var_index = ObType_LookUp(variable)
            
            if var_index == _missing:
                print "pyDart.list:  requested variable:  ",variable," does not exist"
                print "pyDart.list:  Valid variables are:  ", ObType_LookUp(variable,Print_Table=True)
                h5file.close()
                return
            
            if type(variable) == type(1):
                var_index = variable
            
            if self.index == None:
                table.get_where_list("kind == var_index")
            
            if len(self.index) != 0:
                if self.verbose: print "Number of observations:  ",len(self.index)
                data   = self.get_data()
                number = data['number']
                value  = data['value']
                time   = data['date']
                utime  = data['utime']

                lat    = data['lat']
                lon    = data['lon']
                x      = N.where(data['x'] != _missing, data['x']/1000., _missing)
                y      = N.where(data['y'] != _missing, data['y']/1000., _missing)
                z      = N.where(data['z'] != _missing, data['z']/1000., _missing)
                az     = N.where(data['azimuth']   != _missing, data['azimuth'],   _missing)
                el     = N.where(data['elevation'] != _missing, data['elevation'], _missing)
                sat0   = N.where(data['satellite'][:,0] != _missing, data['z']/1000., _missing)
                sat1   = N.where(data['satellite'][:,1] != _missing, data['z']/1000., _missing)
                sat2   = N.where(data['satellite'][:,2] != _missing, data['z']/1000., _missing)

                if dumplength == True:
                    dumplength = len(self.index)
                    print
                    if self.verbose:
                        print("Printing ALL the values, hope it does not take too long because the list has %d entries" % 
                               dumplength)
                    print("%s %s %s" % ("\n", "="*100, "\n"))
                    print "    Index       Variable      Value    Date/Time       Lat    Lon    X(km)  Y(km)  Z(km)      AZ        EL"
                    
                    for n in range(0,dumplength-1):
                        print("%7d     %s     %9.5f    %s  %9.4f  %9.4f  %9.4f  %9.4f  %9.5f  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f" \
                          % (number[n], variable.upper(), value[n], time[n], lat[n], lon[n], x[n], y[n], z[n], az[n], el[n], sat0[n], sat1[n], sat2[n]))
                else:
                    dumplength = min(100,len(self.index))
                    print
                    if self.verbose:
                        print("Printing the first and last 100 values of search indices")
                    print("%s %s %s" % ("\n", "="*100, "\n"))
                    print("    Index     Variable        Value          Date/Time           Lat        Lon        X(km)      Y(km)      Z(km)     AZ      EL")
                    
                    for n in range(0,dumplength-1):
                        print("%7d     %s     %9.5f    %s  %9.4f  %9.4f  %9.4f  %9.4f  %9.5f  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f" \
                          % (number[n], variable.upper(), value[n], time[n], lat[n], lon[n], x[n], y[n], z[n], az[n], el[n], sat0[n], sat1[n], sat2[n]))
                    
                    for n in range(len(time)-dumplength,len(time)):
                        print("%7d     %s     %9.5f    %s  %9.4f  %9.4f  %9.4f  %9.4f  %9.5f  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f" \
                          % (number[n], variable.upper(), value[n], time[n], lat[n], lon[n], x[n], y[n], z[n], az[n], el[n], sat0[n], sat1[n], sat2[n]))
            else:
                print "NO OBSERVATIONS FOUND FOR ", variable.upper()
        
        h5file.close()
        
        return
#-------------------------------------------------------------------------------
    
    def stats(self,variable=None,dumplength=None):
        
        if variable == None:
            if self.verbose: print "\n pyDART.stats:  No variable supplied, not valid, exiting \n"

# Open DART PyTables file
        
        h5file, table = open_pyDart_file(self.hdf5, verbose=self.verbose)
        
                                    # User requests specific variables, search for DART variable index
        var_index = ObType_LookUp(variable)
            
        if var_index == _missing:
            print "pyDart.list:  requested variable:  ",variable," does not exist"
            print "pyDart.list:  Valid variables are:  ", ObType_LookUp(variable,Print_Table=True)
            h5file.close()
            return
            
        if self.index == None:
            table.get_where_list("kind == var_index")
        
        if len(self.index) != 0:
            if self.verbose: print "Number of observations:  ", len(self.index)
            data   = self.get_data()
            number = data['number']
            value  = data['value']
            time   = data['date']
            utime  = data['utime']
            lat    = data['lat']
            lon    = data['lon']
            x      = N.ma.masked_equal(data['x'], 9999.)
            x      = N.ma.masked_equal(data['x'], 9999.)
            y      = N.ma.masked_equal(data['y'], 9999.)
            z      = N.ma.masked_equal(data['z'], _missing)
            z      = N.ma.masked_equal(data['z'], 9999.)
            az     = N.ma.masked_equal(data['azimuth'], _missing)
            el     = N.ma.masked_equal(data['elevation'], _missing)

            x      = N.sort(N.unique(x))
            dx     = N.diff(x)
            
            y      = N.sort(N.unique(y))
            dy     = N.diff(y)
            
            z      = N.sort(z)
            dz     = N.diff(z)
            
            print
            print "%s            Max/Min:   %10.4f   %10.4f  " % (variable, value.max(), value.min())
            print "%s         Mean/stdev:   %10.4f   %10.4f  " % (variable, value.mean(), value.std())
            print "%s  Max/Min    X (m):    %10.1f   %10.1f  " % (variable,x.min(), x.max())
            print "%s  Max/Min   DX (m):    %10.1f   %10.1f  STD: %10.4f " % (variable,dx.max(), dx.min(), dx.std())
            print "%s  Max/Min    Y (m):    %10.1f   %10.1f  " % (variable,y.min(), y.max())
            print "%s  Max/Min   DY (m):    %10.1f   %10.1f  STD: %10.4f " % (variable,dy.max(), dy.min(), dy.std())
            print "%s  Max/Min    Z (m):    %10.1f   %10.1f" % (variable, z.max(), z.min())
            print "%s  Max/Min   Azimuth:   %10.4f   %10.4f" % (variable, az.max(), az.min())
            print "%s  Max/Min Elevation:   %10.4f   %10.4f" % (variable, el.max(), el.min())
            print "%s  Max/Min  Latitude:   %10.4f   %10.4f" % (variable, lat.max(), lat.min())
            print "%s  Max/Min Longitude:   %10.4f   %10.4f" % (variable, lon.max(), lon.min())
            print

        else:
            print "NO OBSERVATIONS FOUND FOR ", variable.upper()
        
        h5file.close()
        
        return

#-------------------------------------------------------------------------------
    
    def ascii2hdf(self,filename=None,radar_loc=None):
        
        if filename != None:
            self.file(filename=filename)
        
        if radar_loc != None:
            print "PyDart.ascii2hdf:  radar location supplied:  lat: %f  lon: %f hgt:  %f" % (radar_loc[0],radar_loc[1], radar_loc[2])

# Create PyTables file
        
        filter_spec = Filters(complevel=5, complib="zlib", shuffle=1, fletcher32=0)
        
        h5file = open_file(self.hdf5, mode = "w", title = version_string, filters=filter_spec)
        
        group_ob_kinds = h5file.create_group("/", 'obs', 'Obs for DART file')

# Create group for observation descriptions
        
        table_ob_kinds = h5file.create_table(group_ob_kinds, 'kinds', DART_ob_kinds, 'Observation Descriptions')
        
        fi = open(self.ascii, 'r')
        fi.readline()                       # Read(str) "obs_sequence"
        fi.readline()                       # Read(str) "obs_kind_definitions"
        ob_kinds = long(fi.readline())      # Read(int) "number of observation types" 
        
        if self.debug:  print 'Number of observation types:  ', ob_kinds
        
        n = 0
        
        row = table_ob_kinds.row
        obtype_dict = {}
        
        while n < ob_kinds:
            stuff = fi.readline()
            stuff = stuff.split()
            row['index'] = int(stuff[0])
            row['name']  = stuff[1]
            print 'Observation kind definitions:  ', row['index'], row['name']
            obtype_dict[int(stuff[0])] = stuff[1]
            row.append()
            n += 1
        
        if self.debug:  print 'Completed reading obs_kind_definitions'

# Create group for misc file header information
        
        group_header = h5file.create_group("/", 'header', 'Header Information for DART file')
        table_header = h5file.create_table(group_header, 'attributes', DART_header, 'Attributes of the observational file')
        
        stuff      = fi.readline()          # Read(str) "num_copies" line
        stuff      = stuff.split()
        
        num_copies = long(stuff[1])         # Important, if you have a truth file, you need to know this number (either 1 or 2)
        num_qc     = long(stuff[3])         
        
        stuff       = fi.readline()         # Read(str) "num_obs" line
        stuff       = stuff.split()
        num_obs     = long(stuff[1])
        max_num_obs = long(stuff[3])
        
        data_storage = []

        for numcp in N.arange(num_copies):
            data_storage.append(fi.readline())  # Read(str) how the data is written....e.g., "observations" or "truth"
            
        for numqc in N.arange(num_qc):          # If the num_qc flag is > 0, read the description of the qc flags
            data_storage.append(fi.readline())

        print "\n ", data_storage
        for numcp in N.arange(num_copies+num_qc):      # Now read the observation, and if provided, the truth value
            print "Observation data types:  %s" % (data_storage[numcp][:-2])
            if ((data_storage[numcp].find("obs") != -1) and 
                (data_storage[numcp].find("tru") != -1) and 
                (data_storage[numcp].find("QC") != -1)
               ):
                print "\n\npyDart.ascii2hdf CANNOT process DART header as observation data storage type is does not have "
                print "%s  %s    %s" % ("'obs'", "'truth'", "'QC'") 
                print "in the description.....EXITING (this can be fixed, see developer)....\n" 
                sys.exit(-1)
        print "\n"
            
        stuff       = fi.readline()         # Read(str) "first   1   last   No of obs" line
        stuff       = stuff.split()
        first       = long(stuff[1])
        last        = long(stuff[3])
        
        row = table_header.row
        
        row['origin_file'] = "Original DART observation file is: " + self.ascii + "\n"
        row['num_copies']  = num_copies
        row['num_qc']      = num_qc
        row['num_obs']     = num_obs
        row['max_num_obs'] = max_num_obs
        row['first']       = first
        row['last']        = last
        
        row.append()
        table_header.flush()
        
        if self.verbose and self.debug:
            print "Number of observation copies:  ", num_copies
            print "Number of QC'd observations:   ", num_qc
            print "Number of observations:        ", num_obs
            print "Max number of observations:    ", max_num_obs

# Find the obs group to create table in
        
        root         = h5file.root
        group_obs    = root.obs
        group_header = root.header

# create table that will hold the observation information
        
        table_obs = h5file.create_table(group_obs, 'observations', DART_obs, 'Observations from DART file')
        
        row = table_obs.row

# Read in the obs sequentially!
        
        n = 0
        
        while True:
            
            stuff            = fi.readline()        # Read(str) "OBS...." line
            if not stuff: break
            
            stuff            = stuff.split()
            
            row["number"]    = long(stuff[1])
            
            for numcp in N.arange(num_copies+num_qc):      # Now read the observation, and if provided, the truth value
                if data_storage[numcp].find("obs") != -1:
                    row["value"]     = read_double_precision_string(fi.readline())
                elif data_storage[numcp].find("tru") != -1:
                    row["truth"]     = read_double_precision_string(fi.readline())
                elif data_storage[numcp].find("QC") != -1:
                    row["qc"]     = read_double_precision_string(fi.readline())
                else:
                    print "pyDART ascii2hdf:  Problem, data_storage is not 'observations' or 'truth' --> exiting!"
                    print data_storage
                    sys.exit(-1)
            
            stuff            = fi.readline()
            stuff            = stuff.replace(","," ").split()
            row["previous"]  = long(stuff[0])
            row["next"]      = long(stuff[1])
            row["cov_group"] = long(stuff[2])
            
            fi.readline()
            fi.readline()
            
            stuff            = fi.readline()
            stuff            = stuff.split()
            row["lon"]       = N.rad2deg(read_double_precision_string(stuff[0]))
            row["lat"]       = N.rad2deg(read_double_precision_string(stuff[1]))
            row["height"]    = read_double_precision_string(stuff[2])

            if row['lon'] > 180.0: row['lon'] = row['lon'] - 360.

            if len(stuff) == 4:
                row["vert_coord"] = long(stuff[3])
            else:
                row["vert_coord"] = long(fi.readline())
            
            fi.readline()
            
            row["kind"]       = int(fi.readline())
            
# Check to see if this GEOS cloud pressure observation
            
            if (row["kind"] == ObType_LookUp("GOES_CWP_PATH") or 
                row["kind"] == ObType_LookUp("GOES_IWP_PATH") or 
                row["kind"] == ObType_LookUp("GOES_LWP_PATH") or 
                row["kind"] == ObType_LookUp("GOES_CWP_ZERO")):
                stuff = fi.readline()
                if stuff.find("2*") > 0: 
                    row["satellite"][0] = stuff.split(" 2*")[1]
                    row["satellite"][1] = stuff.split(" 2*")[1]
                else:
                    stuff = stuff.split(",")
                    row["satellite"][0] = N.rad2deg(read_double_precision_string(stuff[0]))
                    row["satellite"][1] = N.rad2deg(read_double_precision_string(stuff[1]))
                stuff = fi.readline().split()
                row["satellite"][2] = N.float(stuff[0])

# Since pyDart has "standard" integer IDs for observation types, we need to "reset" the "kind" integer
            
            try:
                row["kind"] = ObType_LookUp(obtype_dict[row["kind"]])
            except:
                pass

# Store all geo-coordinates as lat(deg)/lon(deg)/height(m)/az(deg)/el(deg)

# Check to see if its radial velocity and add platform information...need BETTER CHECK HERE!
            
            if row["kind"] == ObType_LookUp("VR"):
                
                fi.readline()
                fi.readline()
                
                stuff                      = fi.readline()
                stuff                      = stuff.split()
                row['platform_lon']        = N.rad2deg(read_double_precision_string(stuff[0]))
                row['platform_lat']        = N.rad2deg(read_double_precision_string(stuff[1]))
                row['platform_height']     = read_double_precision_string(stuff[2])

                if row['platform_lon'] > 180.0: row['platform_lon'] = row['platform_lon'] - 360.

                if len(stuff) == 4:
                    row["platform_vert_coord"] = long(stuff[3])
                else:
                    row["platform_vert_coord"] = long(fi.readline())      
                fi.readline()
                
                stuff                    = fi.readline()
                stuff                    = stuff.replace(","," ")                 
                stuff                    = stuff.split()
                
                row['platform_dir1']     = read_double_precision_string(stuff[0])
                row['platform_dir2']     = read_double_precision_string(stuff[1])
                row['platform_dir3']     = read_double_precision_string(stuff[2])
                row['platform_nyquist']  = read_double_precision_string(fi.readline())
                row['platform_key']      = long(fi.readline())

                row['z']                           = row['height'] - row['platform_height']
                row['x'], row['y'], row['azimuth'] = dll_2_dxy(row['platform_lat'], row['lat'], \
                                                     row['platform_lon'], row['lon'], azimuth=True, degrees=True)

                elevation_angle = N.arcsin(row['platform_dir3']) 

                row['elevation']         = N.rad2deg(elevation_angle)
        
                if elevation_angle > 89.9:   # radar pointing straight up?
                    row['azimuth'] = 0.0

# For reflectivity - lets figure out what the elevation and azimuth are

            if row["kind"] == ObType_LookUp("DBZ"):

# Compute some coordinates 
                if radar_loc != None:
                    radar_lat = radar_loc[0]
                    radar_lon = radar_loc[1]
                    radar_hgt = radar_loc[2]

                    row['z']                           = row['height'] - radar_hgt
                    row['x'], row['y'], row['azimuth'] = dll_2_dxy(radar_lat, row['lat'], radar_lon, row['lon'], azimuth=True, degrees=True)

                    R = N.sqrt(row['x']**2 + row['y']**2)
                    elevation_angle = beam_elv(R, row['z'])
                            
                    row['platform_dir1']     = (row['x'] / R) * N.deg2rad(elevation_angle)
                    row['platform_dir2']     = (row['y'] / R) * N.deg2rad(elevation_angle)
                    row['platform_dir3']     = N.sin(N.deg2rad(elevation_angle))
                    row['platform_nyquist']  = _missing
                    row['platform_key']      = _missing            
                    row['elevation']         = elevation_angle
        
                    if elevation_angle > 89.9:   # radar pointing straight up?
                        row['azimuth'] = 0.0
                
# Okay, now get the time information 

            stuff            = fi.readline()
            stuff            = stuff.split()
            row['seconds']   = int(stuff[0])
            row['days']      = int(stuff[1])

# Add in a full date and time string to data set, as well as create a UTIME in seconds for searching
            
            days = float(stuff[1])+float(stuff[0])/86400.
            date = day_utime.num2date(days)
            
            row['utime']     = round(sec_utime.date2num(date)) #  to prevent sometimes truncating down to next integer
            row['date']      = date.strftime(time_format)
            
            row['error_var'] = read_double_precision_string(fi.readline())
            row['index']     = n
            
            if n % 5000 == 0:
                print "read_DART_ob:  Processed observation # ", n+1, days #,sec_utime.date2num(date) #,py_datetime(start[0],start[1],start[2],start[3],start[4],start[5])
                print "date = ", date
            
            if (N.isnan(row['platform_dir1']) or N.isnan(row['platform_dir2'])):
                print("Found a bad value: %3.3i, skipping\n" % n)
                pass
            else:
                n += 1
                row.append()
            
            table_obs.flush
            
        if n != group_header.attributes.cols.num_obs[0]:
            group_header.attributes.cols.last[0] = n
            group_header.attributes.cols.max_num_obs[0]= n
            group_header.attributes.cols.num_obs[0] = n
            print("Changed number of obs to: %3.3i \n" % (n))
                        
        h5file.close()
        
        fi.close()
        
        print "pyDART.ascii2h5:  Converted ascii DART file to HDF5 DART file"
        
        return 1
#-------------------------------------------------------------------------------
    
    def dict2hdf(self,obs_dict,filename=None,radar_loc=None):
        
        if filename != None:
            self.file(filename=filename)
        
        if radar_loc != None:
            print "PyDart.dict2hdf:  radar location supplied:  lat: %f  lon: %f " %  \
                  (N.degrees(radar_loc[0]),N.degrees(radar_loc[1]))
        
        print
        print "!!! Warning from pyDart.dict2hdf !!!"
        print "pyDART_version_2:  Version 2 now stores lat and lon correctly - please make sure your analysis does as well..."
        print

# Create PyTables file
        
        filter_spec = Filters(complevel=5, complib="zlib", shuffle=1, fletcher32=0)
        
        h5file = open_file(self.hdf5, mode = "w", title = version_string, filters=filter_spec)
        group_ob_kinds = h5file.create_group("/", 'obs', 'Obs for DART file')

# Create group for observation descriptions
        
        table_ob_kinds = h5file.create_table(group_ob_kinds, 'kinds', DART_ob_kinds, 'Observation Descriptions')
        
        row = table_ob_kinds.row
        
        row['index'] = ObType_LookUp("VR")
        row['name']  = "DOPPLER_RADIAL_VELOCITY" #"VR"
        row.append()
        row['index'] = ObType_LookUp("DBZ")
        row['name']  = "RADAR_REFLECTIVITY" #"DBZ"
        row.append()
        
        if self.debug:  print 'Completed reading obs_kind_definitions'

# Create group for misc file header information
        
        group_header = h5file.create_group("/", 'header', 'Header Information for DART file')
        table_header = h5file.create_table(group_header, 'attributes', DART_header, 'Attributes of the observational file')
        
        print table_header
        row = table_header.row
        
        row['origin_file'] = "Original DART observation file is: " + self.ascii + "\n"
        row['num_copies']  = 1
        row['num_qc']      = 1
        row['num_obs']     = obs_dict['counter']
        row['max_num_obs'] = obs_dict['counter']
        row['first']       = 1
        row['last']        = obs_dict['counter']
        
        row.append()
        table_header.flush()
        
        if self.debug:
                print "Number of observation copies:  ", 1
                print "Number of QC'd observations:   ", 1
                print "Number of observations:        ", obs_dict['counter']
                print "Max number of observations:    ", obs_dict['counter']

# Find the obs group to create table in
        
        root         = h5file.root
        group_obs    = root.obs
        group_header = root.header

# create table that will hold the observation information
        
        table_obs = h5file.create_table(group_obs, 'observations', DART_obs, 'Observations from DART file')
        
        row = table_obs.row

# Read in the obs sequentially!
        
        for n in N.arange( obs_dict['counter'] ):
            
            row['number']    = long(n+1)
            row['value']     = obs_dict['value'][n]
            
            if n == 0:
                    row['previous'] = -1
            else:
                    row['previous'] = n-1
            
            row['next']       = n+1
            row['cov_group']  = -1

# Here is where you sometimes need to switch lat and lon in the codes....
            
            row['lon']        = obs_dict['lon'][n]
            row['lat']        = obs_dict['lat'][n]
            row['height']     = obs_dict['height'][n]
            row['vert_coord'] = 3
            row['kind']       = obs_dict['kind'][n]

# Check to see if its radial velocity and add platform information...need BETTER CHECK HERE!
            
            if row['kind'] == 11:

# Need to switch the lat and lon these around to be correct in version 2 of pyDART...
                
                row['platform_lon']        = obs_dict['platform_lon'][n]
                row['platform_lat']        = obs_dict['platform_lat'][n]
                row['platform_height']     = obs_dict['platform_height'][n]
                row['platform_vert_coord'] = obs_dict['platform_vert_coord'][n]
                
                row['platform_dir1']       = obs_dict['platform_dir1'][n]
                row['platform_dir2']       = obs_dict['platform_dir2'][n]
                row['platform_dir3']       = obs_dict['platform_dir3'][n]
                row['platform_nyquist']    = obs_dict['platform_nyquist'][n]
                row['platform_key']        = obs_dict['platform_key'][n]
                
                row['elevation']           = N.rad2deg(N.arcsin(obs_dict['platform_dir3'][n]))
            
                if N.abs(N.cos(row['elevation'])) < 0.001:   # pointing straight up?
                    row['azimuth'] = 0.0
                else:
                    row['azimuth'] = N.rad2deg( N.arctan2( row['platform_dir1'] / N.cos(N.deg2rad(row['elevation'])), \
                                                           row['platform_dir2'] / N.cos(N.deg2rad(row['elevation'])) ) )
                    if row['azimuth'] < 0: row['azimuth'] = row['azimuth'] + 360.
                    
            if row['kind'] == ObType_LookUp("DBZ"):    # See if there is information about the azimuth and elevation supplied
            
                if 'azimuth' in obs_dict.keys():
                    row['azimuth'] = obs_dict['azimuth'][n]
                    if row['azimuth'] < 0: row['azimuth'] = row['azimuth'] + 360.

                if 'elevation' in obs_dict.keys():
                    row['elevation'] = obs_dict['elevation'][n]
                
# Add in information about the observations relative location to radar.  Do this only for radar observation (for now)
            
            if radar_loc != None:
                
                radar_lat = radar_loc[0]
                radar_lon = radar_loc[1]
                                
                row['x'], row['y'] = dll_2_dxy(radar_lat, row['lat'], radar_lon, row['lon'])
                row['z']           = row['height']

# Add in a full date and time string to data set, as well as create a UTIME in seconds for searching
            
            row['seconds']   = obs_dict['sec'][n]
            row['days']      = obs_dict['day'][n]
            
            days = float(obs_dict['day'][n]) + float(obs_dict['sec'][n])/86400.
            date = day_utime.num2date(days)
            
            row['date']      = str(date)
            row['utime']     = round(sec_utime.date2num(date))  # round to prevent sometimes truncating down to next integer
            
            row['error_var'] = obs_dict['error_var'][n]
            row['index']     = n
            
            if n % 5000 == 0:
                print "PyDart.dict2hdf:  Processed observation # ", n+1, days
                print " Date,sec_utime,utime = ", date,sec_utime.date2num(date), row['utime']
            
            row.append()
        
        table_obs.flush
        
        h5file.close()
        
        print "pyDART.dict2h5:  Converted DICTIONARY DATA to HDF5 DART file"
        
        return 1

#-------------------------------------------------------------------------------
# read MRMS created superobbed sweep (on a rectangular grid) and create a
#      DART data base file.  Code also thins zero reflecitivty obs, thresholds
#      a max and min, and removes data close to a radar.  
 
    def mrms(self, ncdf_file, filename=None, lat_bbox=None, lon_bbox=None, hgt_bbox = hgt_bound_box):
    
        if filename != None:
            self.file(filename=filename)

# Open netcdf file(s)

        if ncdf_file.find(".netcdf") != -1:
            fncdf = [ncdf_file]
        else:
            fncdf = glob.glob(ncdf_file+"/*.netcdf")
            
        print "\npyDart/MRMS->  Number of files to be read in:  %d" % (len(fncdf))
        print "\npyDart/MRMS->  First file:  %s " % fncdf[0]
        if len(fncdf) > 1:  print "\npyDart/MRMS->  Last file:  %s " % fncdf[-1]

        if dbz_zeros:
            print("\npyDart/MRMS:  Values of reflectivity < %d are set to clear_air values of: %d" % (dbz_min, dbz_clear_air))
        else:
            print("\npyDart/MRMS:  dbz_zero flag is False")
            
        if dbz_missing_zeros:
            print("\npyDart/MRMS:  Missing reflectivity are set to clear_air values of  %d" % dbz_clear_air)
        else:
            print("\npyDart/MRMS:  dbz_missing zero flag is False")
            
        if dbz_max:
            print("\npyDart/MRMS: Values of reflectivity greater than %d will be set to %d" % (dbz_max, dbz_max))

        print("\npyDart/MRMS:  dBZ values will be thinned horizontally by a factor of %d" % dbz_thin)

        print("\npyDart/MRMS:  dBZ values will be thinned vertically by a factor of %d" % dbz_thin_z)

        if dbz_clear_thin > 1:
             print("\npyDart/MRMS:  Clear air reflectivity values thinned by an additional factor of %d" % dbz_clear_thin)
        else:
             print("\npyDart/MRMS:  No extra thning of clear air data requested")

        print "\npyDart/MRMS->  Lat Bounding Box on:  %f to %f" % (lat_bbox[0], lat_bbox[1])
        print "\npyDart/MRMS->  Lon Bounding Box on:  %f to %f" % (lon_bbox[0], lon_bbox[1])

# Create PyTables file
        
        filter_spec = Filters(complevel=5, complib="zlib", shuffle=1, fletcher32=0)
        
        h5file = open_file(self.hdf5, mode = "w", title = version_string, filters=filter_spec)
        
        group_ob_kinds = h5file.create_group("/", 'obs', 'Obs for DART file')
                
        group_header = h5file.create_group("/", 'header', 'Header Information for DART file')

# Find the obs group to create table in
        
        root         = h5file.root
        group_obs    = root.obs
        group_header = root.header

# create table that will hold the observation information
        
        table_obs = h5file.create_table(group_obs, 'observations', DART_obs, 'Observations from DART file')
    
        row = table_obs.row
        
# set counters

        count       = 0
        count_dbz   = 0
        count_zeros = 0
                  
#  Open netcdf files, read data, flatten, and append

        for file in fncdf:
            
            print "\npyDart/MRMS->  Processing file:  %s" % file

            f = ncdf.Dataset(file, "r")

            lats  = f.variables['Lat'][::dbz_thin]
            lons  = f.variables['Lon'][::dbz_thin]
            
            try:
                hgts  = f.variables['Hgt'][0:-1:dbz_thin_z]
            except KeyError:
                try:
                    hgts  = f.variables['Ht'][0:-1:dbz_thin_z]
                except KeyError:
                    hgts = N.array((dbz_clear_hgts[1],))

  # For radar reflectivity - do some pre-processing based on flags at top of file
  # if dbz_zeros - for all MISSING values, set to zero - except we use the dbz_thin to thin zero reflectivity
  # if dbz_min   - then set all values of dbz < dbz_min to MISSING - this helps create a halo of missing values
  #                around the actual storms.
  # if dbz_max   - finally, as a last step, clip the reflectivity to be less than dbz_max

            dbz = None
           
            try:
                dbz = f.variables['MergedReflectivityQCComposite'][0,::dbz_thin,::dbz_thin]
                print("\n==>MRMS2HDF:  Found MergedReflectivityQC variable in netCDF file!\n")
                dbz = N.expand_dims(dbz, axis=0)
            except KeyError:
                print("\n==>MRMS2HDF:  Cannot find MergedReflectivityQCComposite variable in netCDF file!")

            try:
                dbz = f.variables['MergedReflectivityQC'][0,::dbz_thin_z,::dbz_thin,::dbz_thin]
                print("\n==>MRMS2HDF:  Found MergedReflectivityQC variable in netCDF file!\n")
            except KeyError:
                print("\n==>MRMS2HDF:  Cannot find MergedReflectivityQC variable in netCDF file!")

            try:
                dbz  = f.variables['mrefl_mosaic'][0,::dbz_thin_z,::-dbz_thin,::dbz_thin]
                print("\n==>MRMS2HDF:  Found mrefl_mosaic variable in netCDF file!\n")
            except KeyError:
                print("\n==>MRMS2HDF:  Cannot find mrefl_mosaic variable in netCDF file!\n")

            if dbz == None:
                print("\n==>MRMS2HDF:  Cannot find any of the specified reflectivity variables in netCDF file, exiting....\n")
                raise SystemExit
                         
# Figure out the bounding box indices to make the output go a LOT faster

            index  = (lats > lat_bbox[0] ) & (lats < lat_bbox[1])
            lindex = N.arange(lats.size)[index]
            index  = (lons > lon_bbox[0] ) & (lons < lon_bbox[1])
            mindex = N.arange(lons.size)[index]
            index  = (hgts >= hgt_bbox[0] ) & (hgts < hgt_bbox[1])
            kindex = N.arange(hgts.size)[index]
            
            # print dbz.shape, lats.shape, lons.shape
            # print lindex[0], lats[lindex[0]], lats[lindex[-1]], lindex[-1]
            # print mindex[0], lons[mindex[0]], lons[mindex[-1]], mindex[-1]
            # print kindex, hgts

            for k in kindex:

                if dbz_clear_hgts[0] == True:
                    if hgts[k] == dbz_clear_hgts[1] or hgts[k] == dbz_clear_hgts[2]:
                        hgt_flag = 0
                    else:
                        hgt_flag = 1
                else:
                    hgt_flag = 0

                hgt = hgts[k]       
                        
                for l in lindex:
                    
                    lat = lats[l]
                        
                    for m in mindex:
                        
                        lon = lons[m]
                        
                        thin_zeros = hgt_flag + (l*m) % dbz_clear_thin
                        
                        data = min(dbz[k,l,m], dbz_max)

                        my_mask = False
                        
                        if data == f.MissingData and dbz_missing_zeros and thin_zeros == 0: 
                              data = dbz_clear_air 
                              data_kind = ObType_LookUp("CLEARAIR_REFLECTIVITY")
                              my_mask = True
                              count_zeros += 1

                        if (data > f.MissingData and data < dbz_min) and dbz_zeros and thin_zeros == 0:  
                              data = dbz_clear_air 
                              data_kind = ObType_LookUp("CLEARAIR_REFLECTIVITY")
                              my_mask = True
                              count_zeros += 1
                                                    
                        if data >= dbz_min:  
                              my_mask = True
                              data_kind = ObType_LookUp("REFLECTIVITY")
                              count_dbz += 1

                        if my_mask:   
                            
                            row['number'] = long(count+1)
                            row['value']  = data
            
                            if  count == 0:
                                row['previous'] = -1
                            else:
                                row['previous'] = count-1
            
                            row['next']       = count+1
                            row['cov_group']  = -1

            # Set lat and lon in the codes....
            
                            row['lon']        = lon
                            row['lat']        = lat
                            row['height']     = hgt
                            row['vert_coord'] = 3
                            row['kind']       = data_kind
        
            # Add in a full date and time string to data set, as well as create a UTIME in seconds for searching

                            dt = ncdf.num2date(f.Time, units="seconds since 1970-01-01 00:00:00")
                            gc = ncdf.date2num(dt, units = "days since 1601-01-01 00:00:00")
            
                            row['days']    = N.int(gc)
                            row['seconds'] = (gc - N.int(gc)) * 86400
            
                            row['date']    = dt.strftime(time_format)
                            row['utime']   = f.Time           
            
                            row['error_var'] = dbz_stdev
                    
                            row['index'] = count
                
                            count = count + 1
    
                            if count % 5000 == 0:
                                print "PyDart.MRMS:  Processed observation # %d" %  (count+1)
                                print " Date,sec_utime,utime = ", row['date'], row['utime']
            
                            row.append()
            f.close()
        
        table_obs.flush()
        
        print("\n >%s<" % ("=" * 79))
        print "\n  PyDart.MRMS->  Total number of observations processed observation %d" % long(count)
        print "\n  PyDart.MRMS->  Total number of dBZ > dbz_min: %d" % long(count_dbz)
        print "\n  PyDart.MRMS->  Total number of ZERO dBZs:     %d" % long(count_zeros)
        print("\n >%s<" % ("=" * 79))
        
# Create header information for the table......
        
        table_ob_kinds = h5file.create_table(group_ob_kinds, 'kinds', DART_ob_kinds, 'Observation Descriptions')
        row = table_ob_kinds.row
        
        if count_dbz > 0:
            row['index'] = ObType_LookUp("DBZ")
            row['name']  = "RADAR_REFLECTIVITY" #"DBZ"
            row.append()

        if count_zeros > 0:
            row['index'] = ObType_LookUp("CLEARAIR_REFLECTIVITY")
            row['name']  = "CLEARAIR_REFLECTIVITY" #"DBZ"
            row.append()
              
        table_ob_kinds.flush()

        if self.debug:  print 'Completed reading obs_kind_definitions'

# Create file header information

        table_header = h5file.create_table(group_header, 'attributes', DART_header, 'Attributes of the observational file')
        
        row = table_header.row
                    
        row['num_copies']  = 1
        row['num_qc']      = 0
        row['num_obs']     = count
        row['max_num_obs'] = count
        row['first']       = 1
        row['last']        = count
        
        row.append()
        
        table_header.flush()
        
        if self.debug:
                print "Number of observation copies:  ", 1
                print "Number of QC'd observations:   ", 0
                print "Number of observations:        ", count
                print "Max number of observations:    ", count
        
        h5file.close()
        
        print "\npyDART.MRMS->  Converted MRMS file to HDF5 DART file"
        
        return
#-------------------------------------------------------------------------------
    
    def correct_ens_output(self, ascii=None):
        
        if self.hdf5 == None:
            print "\n pyDart/correct_ens_output:  No HDF5 file name is defined, please add one to the command line... \n"
            return
            
# Open DART PyTables file
        
        h5file, table = open_pyDart_file(self.hdf5, verbose=self.verbose)

# Open ASCII file
        
        if ascii == None:
            fi = open(self.ascii[:-4]+".tmp.out", "w")
        else:
            fi = open(ascii, "w")
        
        if self.index != None:
            fi.write("%d\n" % (len(self.index)) )
            
        else:
            fi.write(" %d\n" % (attr.col('num_obs')[0]))
                            
        if self.debug:  print "pyDart/correct_ens_output:  Completed writing out header information for ascii CORRECT_ENS file"

# If there is no search index defined, then create a temporary one to loop through all rows..
        
        if self.index == None:
            self.index = arange(table.nrows)
        
        n = 0
        for row in table.itersequence(self.index):
            n += 1
                        
            fi.write("    %14.7f   %14.7f    %14.7f   %8.3f \n" % (row["lon"], row["lat"], row["height"], row["value"] ))
            
            if n % 10000 == 0: print "pyDart/hdf2ascii:  Processed observation # ", n
        
        h5file.close()
        
        fi.close()
        
        print "pyDart/correct_ens_output:  Created CORRECT_ENS input file, N = ", n
        
        return 0

#-------------------------------------------------------------------------------
    
    def hdf2ascii(self, ascii=None, obs_error=None):
        
        if self.hdf5 == None:
            print "pyDart/hdf2ascii:  No HDF5 file name is defined, please add one to the command line..."
            return
            
        if obs_error != None:
                print "HEY!!  Changing standard deviation of fields:  ", obs_error, "\n"
                error_dart_fields = []
                for field in obs_error:
                    error_dart_fields.append(ObType_LookUp(field[0],DART_name=True)[0])
        else:
            error_dart_fields = None


# Open DART PyTables file
        
        h5file, table = open_pyDart_file(self.hdf5, verbose=self.verbose)

# Open ASCII file
        
        if ascii == None:
            fi = open(self.ascii[:-4]+".tmp.out", "w")
        else:
            fi = open(ascii, "w")

# Write out header information
        
        fi.write(" obs_sequence\n")
        fi.write("obs_kind_definitions\n")

# Observation types are in the h5file.root.obs.kinds directory
        
        kinds = h5file.root.obs.kinds
        
        fi.write("       %d\n" % N.size(kinds))
        
        for r in kinds.iterrows():
            fi.write("    %d          %s   \n" % (r['index'], r['name']) )
            if self.debug:  print 'pyDart/hdf2ascii:  Written observational types:  ', r
        
        attr = h5file.root.header.attributes

        nobs = attr.col('num_obs')[0]

        fi.write("  num_copies:            %d  num_qc:            %d\n" % (attr.col('num_copies')[0], 1 ))
        
        if self.index != None:
            fi.write(" num_obs:       %d  max_num_obs:       %d\n" % (len(self.index), len(self.index)) )
        else:
            fi.write(" num_obs:       %d  max_num_obs:       %d\n" % (attr.col('num_obs')[0], attr.col('max_num_obs')[0]))
            
        fi.write("observations\n")
        if attr.col('num_copies')[0] == 2:
            fi.write("truth\n")
        fi.write("QC\n")
                
        if self.index != None:
            fi.write("  first:            %d  last:       %d\n" % (1, len(self.index)) )
            if self.debug:
                print "pyDart/hdf2ascii:  Max number of observations:    ", len(self.index)
        else:
            fi.write("  first:            %d  last:       %d\n" % (attr.col('first')[0], attr.col('last')[0]))
            if self.debug:
                print "pyDart/hdf2ascii:  Max number of observations:    ", attr.col('max_num_obs')[0]
        
        if self.debug:
            print "pyDart/hdf2ascii:  Number of observation copies:  ", attr.col('num_copies')[0]
            print "pyDart/hdf2ascii:  Number of QC'd observations:   ", attr.col('num_qc')[0]
        
        if self.debug:  print "pyDart/hdf2ascii:  Completed writing out header information for ascii DART file"

# If there is no search index defined, then create a temporary one to loop through all rows..
        
        if self.index == None:
            self.index = N.arange(table.nrows)
            
#select_table = table.copy(self.index)
        
        n = 0
        for row in table.itersequence(self.index):
            n += 1
            
            fi.write(" OBS            %d\n" % n )
            fi.write("   %20.14f\n" % row["value"]  )
            
            if attr.col('num_copies')[0] == 2:
                fi.write("   %20.14f\n" % row["truth"]  )

            fi.write("   %20.14f\n" % row["qc"] )    # QC flag
            
# Code (from RLT) to output correct index number for each ob in the output DART file

            if n == 1: 
                fi.write(" %d %d %d\n" % (-1, n+1, row["cov_group"]) ) # First obs.
            elif n == len(self.index):
                fi.write(" %d %d %d\n" % (n-1, -1, row["cov_group"]) ) # Last obs.
            else:
                fi.write(" %d %d %d\n" % (n-1, n+1, row["cov_group"]) ) 
            
            fi.write("obdef\n")
            fi.write("loc3d\n")
 
            if row["lon"] < 0.:
                ob_lon = row["lon"] + 360.
            else:
                ob_lon = row["lon"] 
            
            fi.write("    %20.14f          %20.14f          %20.14f     %d\n" % 
                    (N.radians(ob_lon), N.radians(row["lat"]), row["height"], row["vert_coord"]) )
            
            fi.write("kind\n")
            
            fi.write("     %d     \n" % row["kind"] )

# If this GEOS cloud pressure observation, write out extra information (NOTE - NOT TESTED FOR HDF2ASCII LJW 04/13/15)

            if (row["kind"] == ObType_LookUp("GOES_CWP_PATH") or 
                row["kind"] == ObType_LookUp("GOES_IWP_PATH") or 
                row["kind"] == ObType_LookUp("GOES_LWP_PATH") or 
                row["kind"] == ObType_LookUp("GOES_CWP_ZERO")):

                fi.write("    %20.14f          %20.14f  \n" % (row["satellite"][0], row["satellite"][1]) )
                fi.write("    %20.14f  \n" % (row["satellite"][2]) )

# Check to see if its radial velocity and add platform information...need BETTER CHECK HERE!
            
            if row["kind"] == ObType_LookUp("VR"):
                
                fi.write("platform\n")
                fi.write("loc3d\n")

                if row["platform_lon"] < 0.:
                    plat_lon = row["platform_lon"] + 360.
                else:
                    plat_lon = row["platform_lon"] 

                platform_lat, platform_lon = N.radians(row["platform_lat"]), N.radians(plat_lon)
                
                fi.write("    %20.14f          %20.14f        %20.14f    %d\n" %
                        (platform_lon, platform_lat, row["platform_height"], row["platform_vert_coord"]) )
                
                fi.write("dir3d\n")
                
                dir1, dir2, dir3 = (row["platform_dir1"], row["platform_dir2"], row["platform_dir3"])
                fi.write("    %20.14f          %20.14f        %20.14f\n" % (dir1, dir2, dir3) )
                fi.write("     %20.14f     \n" % row["platform_nyquist"] )
                fi.write("     %d     \n" % row["platform_key"] )

# Done with special radial velocity obs back to dumping out time, day, error variance info
            
            fi.write("    %d          %d     \n" % (row["seconds"], row["days"]) )

# Logic for command line override of observational error variances

            if error_dart_fields != None:
            
                try:
                    std_dev = float(obs_error[error_dart_fields.index(row['kind'])][1])
                    variance = std_dev*std_dev
                    fi.write("    %20.14f  \n" % variance )
                except ValueError:
                    fi.write("    %20.14f  \n" % row['error_var'] )
                
            else:
                fi.write("    %20.14f  \n" % row['error_var'] )

            if n % 10000 == 0: print "pyDart/hdf2ascii:  Processed observation # ", n
        
        h5file.close()
        
        fi.close()
        
        print "pyDart/hdf2ascii:  Created ascii DART file, N = ", n
        
        return 0

#-------------------------------------------------------------------------------
    
    def getDartTimes(self, output_file_name=None):

# Open DART PyTables file
        
        h5file, table = open_pyDart_file(self.hdf5, verbose = self.verbose)
        
        n = 0
        date = []
        nobs = []
        for row in table.iterrows():
            if n == 0:
                date.append(row['date'])
                n = 1
            else:
                if date[n-1] != row['date']:
                    date.append(row['date'])
                    n += 1
        
        if output_file_name == None:
            print n
            for d in date:
                print "%s  %s  %s  %s  %s  %s" % (d[0:4],d[5:7],d[8:10],d[11:13],d[14:16],d[17:19])
        else:
            fi = open(output_file_name, "w")
            fi.write("   %d\n" % n )
            for d in date:
                fi.write("%s  %s  %s  %s  %s  %s\n" % (d[0:4],d[5:7],d[8:10],d[11:13],d[14:16],d[17:19]) )
            fi.close()
        
        h5file.close()
        
        return 0

#-------------------------------------------------------------------------------

def read_double_precision_string(the_string):  # We have a double precision scientific string which I can only read this way

      if the_string.count("D") > 0:
#           match      = re.search('-?\d*(\.\d+|)[Dd][+\-]\d\d\d?', the_string)
#           number_str = str(match.group(0))
#       number     = number_str.replace('D','e')
            number     = the_string.replace('D','e')
      else:
            number     = the_string

      return N.double(number)

#===============================================================================        

# DART FILE HEADER FORMAT  (slashes mean file text...)
#-------------------------------------------------------------------------------
#  /obs_sequence/
#  /obs_kind_defnitions/
#  ob_kinds[int8]
#  ob_kind_definition(1)[int8,str]
#  ob_kind_definition(2)[int8,str]
#       .
#       .
#       .
#       .
#  ob_kind_definition(ob_kinds)[int8,str]
#  /num_copies:/  num_copies[int8]   /num_qc:/  num_qc[int8]
#  /num_obs:   /  num_obs[int8]      /max_num_obs:/  max_num_obs[int8]
#  /observations/
#  /first:/   first[int8]    /last:/  last[int8]
#
#-------------------------------------------------------------------------------
# DART FILE OBSERVATION FORMAT  (slashes mean file text...)
#-------------------------------------------------------------------------------
#   /OBS/   ob_number[int8]
#   value[flt64]
#   prev[int8]  next[int8]  covariance_group[int8]
#   /obdef/
#   /loc3d/
#   ob_lat[flt64]  ob_lon[flt64] ob_height[flt64] ob_vert_coord[int8]
#   /type/
#   ob_type[int8]
#   /platform/              {if ob_type == Doppler_velocity, need to know where radar is...)
#   /loc3d/
#   plat_lat[flt64]  plat_lon[flt64] plat_height[flt64] plat_vert_coord[int8]
#   /dir3d/
#   dir1[flt64], dir2[flt64], dir3[flt64]
#   nyquist[flt32]          {Dowell convention adds in Nyquist information}
#   ob_def_key[int8]
#   seconds[int8]  days[int8]
#   ob_error_variance[flt32]
#

class DART_ob_kinds(IsDescription):
        index               = Int64Col(dflt=long(_missing))
        name                = StringCol(255)

class DART_header(IsDescription):
        origin_file         = StringCol(255)
        num_copies          = Int64Col(dflt=long(_missing))
        num_qc              = Int64Col(dflt=long(_missing))
        num_obs             = Int64Col(dflt=long(_missing))
        max_num_obs         = Int64Col(dflt=long(_missing))
        first               = Int64Col(dflt=long(_missing))
        last                = Int64Col(dflt=long(_missing))

class DART_obs(IsDescription):
        number              = Float64Col(dflt=_missing)
        value               = Float64Col(dflt=_missing)
        truth               = Float64Col(dflt=_missing)
        previous            = Float64Col(dflt=_missing)
        qc                  = Float64Col(dflt=_missing)
        next                = Float64Col(dflt=_missing)
        cov_group           = Float64Col(dflt=_missing)
        lat                 = Float64Col(dflt=_missing)
        lon                 = Float64Col(dflt=_missing)
        height              = Float64Col(dflt=_missing)
        vert_coord          = Float64Col(dflt=_missing)
        kind                = Int64Col  (dflt=long(_missing))
        elevation           = Float64Col(dflt=_missing)
        azimuth             = Float64Col(dflt=_missing)
        error_var           = Float64Col(dflt=_missing)
        seconds             = Int64Col  (dflt=long(_missing))
        days                = Int64Col  (dflt=long(_missing))
        date                = StringCol (25)
        utime               = Int64Col  (dflt=long(_missing))
        index               = Int64Col  (dflt=long(_missing))
        x                   = Float64Col(dflt=_missing)
        x_m                 = Float64Col(dflt=_missing)
        y                   = Float64Col(dflt=_missing)
        y_m                 = Float64Col(dflt=_missing)
        z                   = Float64Col(dflt=_missing)
        d                   = Float64Col(dflt=_missing)
        platform_lat        = Float64Col(dflt=_missing)
        platform_lon        = Float64Col(dflt=_missing)
        platform_height     = Float64Col(dflt=_missing)
        platform_vert_coord = Int64Col  (dflt=long(_missing))
        platform_dir1       = Float64Col(dflt=_missing)
        platform_dir2       = Float64Col(dflt=_missing)
        platform_dir3       = Float64Col(dflt=_missing)
        platform_nyquist    = Float64Col(dflt=_missing)
        platform_key        = Int64Col  (dflt=long(_missing))
        Hxb_bar             = Float64Col(dflt=_missing)
        Hxbm                = Float64Col(shape=(30,), dflt=_missing)
        sdHxb               = Float64Col(dflt=_missing)
        Hxa_bar             = Float64Col(dflt=_missing)
        Hxam                = Float64Col(shape=(30,), dflt=_missing)
        sdHxa               = Float64Col(dflt=_missing)
        dep                 = Float64Col(dflt=_missing)
        Yb_prime            = Float64Col(shape=(30,), dflt=_missing)
        satellite           = Float64Col(shape=(3,), dflt=_missing)

#-------------------------------------------------------------------------------
# Main function defined to return correct sys.exit() calls

def main(argv=None):
    if argv is None:
           argv = sys.argv

# Init some local variables
    
    start = None
    end   = None

# Initialize class object...
    
    myDART = pyDART(verbose=_verbose, debug=_debug)

# Command line interface for PyDart
    
    parser = OptionParser()
    parser.add_option("-f", "--file",        dest="file",      type="string", help = "Filename of ascii or PyDART file for conversion/search")
    parser.add_option("-l", "--list",        dest="list",      default=False, help = "Boolean flag to list basic contents of the file",            action="store_true")
    parser.add_option(      "--stats",       dest="stats",     default=False, help = "Gives basic stats for variable, helpful to compare files",   action="store_true")
    parser.add_option("-d", "--dir",         dest="dir",       default=None,  nargs=2, type="string", help = "Directory of files to process and file suffix [*.out, *VR.h5]") 
    parser.add_option(      "--ascii2hdf",   dest="ascii2hdf", default=False, help = "Boolean flag to convert ascii DART file to HDF5 DARTfile",   action="store_true")

    parser.add_option(      "--hdf2ascii",   dest="hdf2ascii", default=False, help = "Boolean flag to convert HDF5 DART file to ascii DART file",  action="store_true")
    parser.add_option(      "--nc2hdf",      dest="nc2hdf",    type="string", help = "File name of file or directory to convert netcdf W2 files to HDF5-DART")
    parser.add_option(      "--mrms",        dest="mrms",      type="string", help = "File name of file or directory to convert netCDF MRMS files to HDF5-DART")
    parser.add_option(      "--start",       dest="start",     type="string", help = "Start time of search in YYYY,MM,DD,HH,MM,SS")
    parser.add_option(      "--end",         dest="end",       type="string", help = "End time of search in YYYY,MM,DD,HH,MM,SS")
    parser.add_option(      "--getDartTimes",dest="DartTimes", default=False, help = "Boolean flag to dump out observations times as in getDARTtimes", action="store_true")
    parser.add_option(      "--condition",   dest="condition", default=None,  type = "string", help = "string having following syntax for searches:  '( z1 < height < z2 )'" )
    parser.add_option(      "--variable",    dest="variable",  default=None,  type = "string", help = "String containing the type of observation to list information:  VR, DBZ")
    parser.add_option(      "--plot",        dest="plot",      default=False, help = "Boolean flag to plot observations data", action="store_true")
    parser.add_option(      "--quiet",       dest="verbose",   default=True,  help = "Boolean flag to suppress text print (default=False)", action="store_false")
    parser.add_option(      "--nodebug",     dest="nodebug",   default=False, help = "Boolean flag to suppress debug output (default=False)", action="store_true")
    parser.add_option("-o", "--output",      dest="output",    default=None,  type = "string", help = "Output filename to store output from various programs")
    parser.add_option(      "--radar_loc",   dest="radar_loc", default=None,  type = "float",  nargs=3, help = "The lat/lon/hgt radar location used in creating X/Y offsets.  Usage:  --radar_loc 35.634 -95.787 370.0")
    parser.add_option("-s", "--search",      dest="search",    default=False, help = "Boolean flag to force search of table - use when search does not use time indices", action="store_true")
    parser.add_option(      "--xloc",        dest="xloc",      default=None,  type = "float",  nargs=2, help = "Search for obs within these x limits. Usage:  --xloc Xmin Xmax (in km)")
    parser.add_option(      "--yloc",        dest="yloc",      default=None,  type = "float",  nargs=2, help = "Search for obs within these y limits. Usage:  --yloc Ymin Ymax (in km)")
    parser.add_option(      "--zloc",        dest="zloc",      default=None,  type = "float",  nargs=2, help = "search for obs within these z limits. Usage:  --zloc Zmin Zmax (in km)")
    parser.add_option(      "--lat_box",     dest="lat_box",   default=None,  type = "float",  nargs=2, help = "Search for MRMS within these lat limits. Usage:  --lat_box lat_south lat_north")
    parser.add_option(      "--lon_box",     dest="lon_box",   default=None,  type = "float",  nargs=2, help = "Search for MRMS within these lon limits. Usage:  --lon_box lon_west lon_east")
    parser.add_option(      "--obserror",    dest="obserror",  default=None,  type = "string", nargs=2, action="append", help = "Change the stored standard deviation of a observational type. Usage: --obserror DBZ 3.0")   
    parser.add_option(      "--merge",       dest="merge",     default=False, help = "Boolean flag to merge several HDF5 obs_seq files", action="store_true")
    parser.add_option(      "--sort",       dest="sort",     default=False, help = "Boolean flag to sort pyDart table in ascending order", action="store_true")
    parser.add_option(      "--correctens",  dest="correctens",default=False, help = "Boolean flag to dump out observed reflectivity to be ingested into correct_ensemble", action="store_true")   
    parser.add_option(      "--addindex",    dest="addindex",  default=False, help = "Boolean flag to create indices for faster search", action="store_true")   
    parser.add_option(      "--scatter",     dest="scatter",   default=False, help = "Boolean flag to scatterplot observations data", action="store_true")
    parser.add_option(      "--rad2deg",     dest="rad2deg",   default=False, help = "Boolean flag to convert old format lat/lon in radians to degrees for plotting", action="store_true")
    
    (options, args) = parser.parse_args()
    
    if options.radar_loc != None:
        radar_loc = [options.radar_loc[0],options.radar_loc[1], options.radar_loc[2]]
    else:
        radar_loc = _radar_loc
    
    if options.rad2deg == True:
        _convertLatLon = True
    else:
        _convertLatLon = False

    if options.nodebug:
        print "Turn off debug"
        myDART.debug = False
    
    if options.verbose:
        myDART.verbose = options.verbose
        print '------------------------------------------------------------------------------------'
        print
        print '  ==> BEGIN PROGRAM PYDART ',version_string
        print
        print '------------------------------------------------------------------------------------'
    
    in_filenames  = []

    if options.dir == None:

        if options.file == None:
            print "\n                NO INPUT FILE or DIRECTORY IS SUPPLIED, EXITING.... \n "
            parser.print_help()
            print
            sys.exit(1)

        else:
            in_filenames.append(os.path.abspath(options.file))

    else:

        suffix = options.dir[1]
        in_filenames = glob.glob("%s/%s" % (os.path.abspath(options.dir[0]), suffix))
        print("\n pyDart:  Processing %d files in the directory:  %s\n" % (len(in_filenames), options.dir[0]))
        print("\n pyDart:  First file is %s\n" % (in_filenames[0]))
        print("\n pyDart:  Last  file is %s\n" % (in_filenames[-1]))

    if options.ascii2hdf:

        for file in in_filenames:
            if file[-3:] == "out":
                if options.verbose:
                    print "\n  PyDart:  converting ASCII DART file:  ", file
                myDART.file(filename = file)
                myDART.ascii2hdf(radar_loc = radar_loc)
                print "\n PyDart:  Completed convertion, PyDART file:  ", myDART.hdf5
            else:
                print "\n File is not labeled `.out`, please rename file..."
    
    if options.start != None:           # Convert start flag to tuple for search
        char = options.start.split(",")
        if len(char) < 5:
            print "pyDart:  Incorrect date/time format for search start (need YYYY,MM,DD,HH,MM,SS)"
            print "pyDart:  submitted arg = ", char, " Exiting program"
            sys.exit(1)
        elif len(char) == 5:
           start =  (int(char[0]), int(char[1]), int(char[2]), int(char[3]), int(char[4]),00)
        else:
           start =  (int(char[0]), int(char[1]), int(char[2]), int(char[3]), int(char[4]), int(char[5]))
        if options.end == None:
            end = (2040,01,01,00,00,00)
        options.search == True
    
    if options.end != None:              # Convert end flag to tuple for search
        char = options.end.split(",")
        if len(char) < 5:
            print "pyDart:  Incorrect date/time format for search end (need YYYY,MM,DD,HH,MM,SS)"
            print "pyDart:  submitted arg = ", char, " Exiting program"
            sys.exit(1)
        elif len(char) == 5:
           end =  (int(char[0]), int(char[1]), int(char[2]), int(char[3]), int(char[4]),00)
        else:
           end =  (int(char[0]), int(char[1]), int(char[2]), int(char[3]), int(char[4]), int(char[5]))
        if options.start == None:
            start = (1970,01,01,00,00,00)
        options.search == True
    
    if options.start == None and options.end == None and not options.ascii2hdf:              # Convert end flag to tuple for search
        start = (1970,01,01,00,00,00)
        end   = (2040,01,01,00,00,00)
        options.search == True
            
    loc =[]
    if options.xloc != None:
        xloc = list(options.xloc)        # convert to list
        xloc = [x * 1000. for x in xloc] # convert to meters (ugly, but the only for lists)
        xloc.sort()                      # make sure the list is from xmin --> xmax...
        loc.append( "(" + str(xloc[0]) + " <= x )" )
        loc.append( "( x <= " + str(xloc[1]) + ")" )
        options.search = True
    
    if options.yloc != None:
        yloc = list(options.yloc)        # convert to list
        yloc = [y * 1000. for y in yloc] # convert to meters (ugly, but the only for lists)
        yloc.sort()                      # make sure the list is from ymin --> ymax...
        loc.append( "(" + str(yloc[0]) + " <= y )" )
        loc.append( "( y <= " + str(yloc[1]) + ")" )
        options.search = True
        
    if options.zloc != None:
        zloc = list(options.zloc)        # convert to list
        zloc = [z * 1000. for z in zloc] # convert to meters (ugly, but only way for lists)
        zloc.sort()                      # make sure the list is from zmin --> zmax...
        loc.append( "(" + str(zloc[0]) + " <= z )" )
        loc.append( "( z <= " + str(zloc[1]) + ")" )
        options.search = True
        
    if options.sort and not options.merge:   # Do a search and return the index of that search
        if options.file[-2:] == "h5":
            sortTable(options.file)
            if options.verbose:
                print("\n pyDart: %s is now sorted by time\n" % options.file)
        else:
            print("\n  pyDart:  ERROR!!  Can only sort on HDF5 pyDART file, exiting...")
            sys.exit(1)
            
    if options.search:   # Do a search and return the index of that search
        if options.file[-2:] == "h5":
            myDART.file(filename = options.file)
            myDART.search(variable=options.variable, start = start, end = end, condition=options.condition, loc=loc)
            if options.verbose:
                if myDART.index != None:
                    print("\n pyDart: %d Observations found between %s and %s in file %s " % \
                          (len(myDART.index), str(start), str(end), myDART.hdf5))
                else:
                    print("\n pyDart:  No Observations found between %s and %s in file %s " % (str(start), str(end), myDART.hdf5))
                    sys.exit(1)
        else:
            print("\n  pyDart:  ERROR!!  Can only search on HDF5 pyDART file, exiting...")
            sys.exit(1)
    
    if options.list and options.variable == None:
        if myDART.verbose:  print("\n PyDart:  Listing file contents")
        myDART.file(filename = options.file)
        myDART.list()

    if options.list and options.variable != None:
        if myDART.verbose:  print("\n PyDart:  Listing information requested about variable:  ", options.variable)
        myDART.file(filename = options.file)
        myDART.list(variable = options.variable, dumplength = True)
    
    if options.stats:
        if myDART.verbose:  print("\n PyDart:  Creating stats")
        myDART.file(filename = options.file)
        myDART.stats(variable = options.variable)

    if options.plot:
        if myDART.verbose: print("\n PyDart:  plotting data")
        myDART.file(filename = options.file)
        myDART.plot(variable = options.variable, convertLatLon = _convertLatLon, savefig=options.output)

    if options.scatter:
      print("\n PyDart:  making scatterplot of data")
      myDART.file(filename = options.file)
      myDART.scatter(variable = options.variable, convertLatLon = _convertLatLon)
    
    if options.hdf2ascii:

        for file in in_filenames:
            if myDART.index != None:
                if myDART.verbose:  print("\n PyDart:  converting HDF5 DART file:  %s" % file)
                myDART.file(filename = file)
                myDART.hdf2ascii(obs_error=options.obserror)
                if myDART.verbose:  print("\n PyDart:  Completed convertion, PyDART file:  %s" % myDART.ascii)
            else:
                print("\n PyDart:  No search indices supplied, so converting entire h5 file to ascii...")
                myDART.file(filename = file)
                myDART.hdf2ascii(obs_error=options.obserror)

    if options.nc2hdf:
        myDART.file(filename = options.file)
        myDART.nc2hdf(options.nc2hdf)
        if myDART.verbose:  print("\n PyDart:  Completed convertion, PyDART file:  %s" % options.ncdf)

    if options.mrms:
        if options.lat_box:  
            lat_bbox = options.lat_box
        else:
            lat_bbox = lat_bound_box

        if options.lon_box:  
            lon_bbox = options.lon_box
        else:
            lon_bbox = lon_bound_box
            
        myDART.mrms(options.mrms, filename = options.file, lat_bbox = lat_bbox, lon_bbox = lon_bbox )
        myDART.hdf2ascii()
        
        if myDART.verbose:  print("\n PyDart:  Completed convertion, PyDART file:  %s" % options.file)
 
    if options.addindex:
        myDART.file(filename = options.file)
        myDART.indexrows
        if myDART.verbose:  print("\n PyDart:  Completed convertion, PyDART file:  %s" % options.file)
       
    if options.merge and options.file:
        if len(in_filenames) > 2:
            mergeTables(options.file, in_filenames)
            if options.sort:
                sortTable(options.file)
                if options.verbose:
                    print("\n pyDart: %s is now sorted by time" % options.file)
        else:
            print("\n You need to specify a directory and a grep pattern for the two files you want merged using \
                   '-d dir grep_pattern' \n")
            print("\n You need to specify a directory and a grep pattern for the two files you want merged using \
                  '-d dir grep_pattern' \n")
            sys.exit(-1)
   
    if options.correctens:
        if options.verbose and myDART.verbose:  print("\n PyDart:  Creating correct_ensemble input file from observations")
        if options.file[-2:] == "h5":
            myDART.file(filename = options.file)
            myDART.search(variable="DBZ", start = start, end = end, condition=options.condition, loc=loc)
            if options.verbose:
                if myDART.index != None:
                    print("\n pyDart: %d Observations found between %s and %s in file %s " % \
                         (len(myDART.index), str(start), str(end), myDART.hdf5))
                else:
                    print("\n pyDart:  No Observations found between %s and %s in file %s " % (str(start), str(end), myDART.hdf5))
                    sys.exit(0)
            myDART.correct_ens_output()
            if options.verbose and myDART.verbose:  
                print("\n PyDart:  Created correct_ensemble input file from dbz observations \n")
        else:
            print("\n pyDart:  ERROR!!  Can only search on HDF5 pyDART file, exiting...")
            sys.exit(-1)

#-------------------------------------------------------------------------------
# Main program for testing...
#
if __name__ == "__main__":
    sys.exit(main())

# End of file
