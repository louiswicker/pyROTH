#############################################################
# pyROTH: A program to process LVL-2 volumes, unfold radial #
#       velocities, and thin the radardata using a          #
#       Cressman analysis.  Gridded data is then outputted  #
#       to DART format for assimilationp.                   #
#                                                           #
#       Python package requirements:                        #
#       ----------------------------                        #
#          Numpy                                            #
#       Scipy                                               #
#       matplotlib                                          #
#       ARM-DOE Py-ART toolkit                              #
#       PyResample toolkit                                  #
#                                                           #
# Originally coded by Blake Allen, August 2016              #
#############################################################
#
# Big modifications by Lou Wicker September 2016
#
#############################################################
import os
import sys
import glob
import time as timeit

import numpy as np
import scipy.interpolate
import scipy.ndimage as ndimage
import scipy.spatial

import netCDF4 as ncdf

from pyresample import kd_tree, utils, geometry
import pyart
from optparse import OptionParser

import pylab as P  
from mpl_toolkits.basemap import Basemap
from pyart.graph import cm

# Parameter dict for Gridding
_grid_dict = {
              'grid_spacing_xy' : 6000.,   # meters
              'grid_radius_xy'  : 150000., # meters
              'weight_func'     : 'Cressman',
              'ROI'             : 3000./0.707, # meters
              'projection'      : 'lcc', # map projection to use for gridded data
              'mask_vr_with_dbz': True,
              '0dbz_obtype'     : True,
              'halo_footprint'  : 4,
              'nthreads'        : 1,
             }

# Dict for the standard deviation of obs_error for reflectivity or velocity (these values are squared when written to DART)            
_obs_errors = {
                'reflectivity': 5.0,
                'velocity': 2.0
              }

# Parameter dict setting radar data parameters
               
_radar_parameters = {
                     'min_dbz_analysis': 10.0, 
                     'max_range': 150000.,
                    }
        
# This flag adds an k/j/i index to the DART file, which is the index locations of the gridded data
_write_grid_indices = True

# True here uses the basemap county database to plot the county outlines.
_plot_counties = True

#=========================================================================================
# Class variable used as container

class Gridded_Field(object):
  
  def __init__(self, name, data=None, **kwargs):    
    self.name = name    
    self.data = data    
    
    if kwargs != None:
      for key in kwargs:  setattr(self, key, kwargs[key])
      
  def keys(self):
    return self.__dict__

#=========================================================================================
# DBZ Mask

def dbz_masking(ref, thin_zeros=2):

  mask = (ref.data.mask == True)
  
  ref.data.mask = False
    
  ref.data[mask] = 0.0
  
  nlevel = ref.data.shape[0]
  
  for n in np.arange(nlevel):
  
    max_values = ndimage.maximum_filter(ref.data[n,:,:], size=_grid_dict['halo_footprint'])
    halo  = (max_values > 0.1) & mask[n]
    
    ref.data.mask[n,halo] = True

  if thin_zeros > 0:

    for n in np.arange(nlevel):
    
      mask1 = np.logical_and(ref.data[n] < 0.1, ref.data.mask[n] == False)  # true for dbz=0
      mask2 = ref.data.mask[n]            # true for dbz>10
      mask1[::thin_zeros, ::thin_zeros] = False
      ref.data.mask[n] = np.logical_or(mask1, mask2)
    
  return ref
    
#=========================================================================================
# DART obs definitions (handy for writing out DART files)

def ObType_LookUp(name,DART_name=False,Print_Table=False):
      """ObType_LookUp returns the DART kind number for an input variable type.  There seems
         to be several ways observations names are written in the observation inputs
         and output files, e.g., in the DART ascii files and in the ***.obs.nc files,
         so this function is designed to handle the variety of cases and return the
         integer corresponding to the DART definitionp.
      
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
                      "UNFOLDED VELOCITY":                [11,   "DOPPLER_RADIAL_VELOCITY"] ,
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
                      "REFL":                             [12,   "REFLECTIVITY"],
                      "FLASH_RATE_2D":                    [2014, "FLASH_RATE_2D"],
                      "METAR_ALTIMETER":                  [71,   "METAR_ALTIMETER"],
                      "LAND_SFC_ALTIMETER":               [70,   "LAND_SFC_ALTIMETER"],
                      "LAND_SFC_DEWPOINT":                [58,   "LAND_SFC_DEWPOINT"],
                      "LAND_SFC_TEMPERATURE":             [25,   "LAND_SFC_TEMPERATURE"],
                      "LAND_SFC_U_WIND_COMPONENT":        [23,   "LAND_SFC_U_WIND_COMPONENT"],
                      "LAND_SFC_V_WIND_COMPONENT":        [24,   "LAND_SFC_V_WIND_COMPONENT"],
                      "GOES_CWP_PATH":                    [80,   "GOES_CWP_PATH"]
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

########################################################################

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

        elvrad = np.arctan((hgtdb*np.cos(rngdb) - frthrde)/(hgtdb * np.sin(rngdb)))

        return np.rad2deg(elvrad)

    else:

        return -999.

####################################################################################### 
#
# write_DART_ascii is a program to dump radar data to DART ascii files.
#
# Usage: 
#
#   obs:  a gridded data object (described below)
#
#   fsuffix:  a string containing a label for the radar fields
#             the string "obs_seq" will be prepended, and
#             the filename will have ".txt" appended after fsuffix.
#             If fsuffix is not supplied, data will be write to
#             "obs_seq.txt".
#
#   obs_error:  the observation error for the field (2 m/s, 5 dbz)
#               this is NOT the variance, the stddev!
#               YOU MUST SPECIFY the obs_error, or program will quit.
#
#   obs object spec:  The obs object must have the following  attributes...
#
#       obs.data:       3D masked numpy array of radar data on a grid.
#                       missing data are filled with  np.nan, and mask=True
#                       for those points.
#
#       obs.data.mask:  the 3D mask for data array (part of np.ma.array)
#
#       obs.field:  string name of field (valid:  "reflectivity", "velocity") 
#
#       obs.lats/lons:  2D numpy arrays of lats and lons of horizontal grid 
#                       locations in degrees.
#
#       obs.hgts:  3D array of heights in meters AGL.
#
#       obs.radar_hgt:  height of radar in meters above MSL
#
#       obs.time:  a dictorary having {'data': an array containing time in seconds, size could be=1  
#                                      'units': a calender reference for datatime calcs}
#                  example:  data: array([   0.535,    0.56 ,    0.582,...])
#                            units: 'seconds since 2013-05-31T23:55:56Z'
# 
#  -->  for radial velocity, more metadata is needed
#
#        obs.nyquist:  a 1D array of dimension z, containing nyquist velo for each tilt
#
#        obs.radar_lat:  latitude of radar
#        obs.radar_lon:  longitude of radar
#
#        obs.xg:  1D array of x-distance of grid point from radar
#        obs.yg:  1D array of y-distance of grid point from radar
#
#
########################################################################################  
def write_DART_ascii(obs, fsuffix=None, obs_error=None, zero_dbz_obtype=True):

  if fsuffix == None:
      print("\n write_DART_ascii:  No output file name is given, writing to %s" % "obs_seq.txt")
      filename = "obs_seq.out"
  else:
      filename = "%s_%s.out" % ("obs_seq", fsuffix)
      
  if obs_error == None:
      print "write_DART_ascii:  No obs error defined for observation, exiting"
      raise SystemExit

# Open ASCII file for DART obs to be written into.  We will add header info afterward
  
  fi = open(filename, "w")
  
# This is why I love python - they think of everything.  Here Numpy has an interator over 
#      an array, and it will extract the indices for you automatically..so create a single 
#      loop over a MD array is very simple....

  print("\n Writing %s to file...." % obs.field.upper())
    
  data       = obs.data
  lats       = np.radians(obs.lats)
  lons       = np.radians(obs.lons)
  hgts       = obs.zg + obs.radar_hgt
  vert_coord = 3
  kind       = ObType_LookUp(obs.field.upper())
  truth      = 1.0  # dummy variable

# platform information

  if kind == ObType_LookUp("VR"):
      platform_nyquist    = obs.nyquist
      platform_lat        = np.radians(obs.radar_lat)
      platform_lon        = np.radians(obs.radar_lon)
      platform_hgt        = obs.radar_hgt
      platform_key        = 1
      platform_vert_coord = 3

# Use the volume mean time for the time of the volume
      
  dtime   = ncdf.num2date(obs.time['data'].mean(), obs.time['units'])
  days    = ncdf.date2num(dtime, units = "days since 1601-01-01 00:00:00")
  seconds = np.int(86400.*(days - np.floor(days)))
  
# Print the number of value gates

  print("\n Number of good observations:  %d" % np.sum(data.mask[:]==False))
  
# Create a multidimension iterator to move through 3D array creating obs
 
  it = np.nditer(data, flags=['multi_index'])
  nobs = 0
  nobs_clearair = 0
  
  while not it.finished:
      k = it.multi_index[0]
      j = it.multi_index[1]
      i = it.multi_index[2]
      
      if data.mask[k,j,i] == True:   # bad values
          pass
      else:          
          nobs += 1
  
          if _write_grid_indices:
              fi.write(" OBS            %d     %d     %d    %d\n" % (nobs,k,j,i) )
          else:
              fi.write(" OBS            %d     %d     %d    %d\n" % (nobs,k,j,i) )
              
          fi.write("   %20.14f\n" % data[k,j,i]  )
          fi.write("   %20.14f\n" % truth )
            
          if nobs == 1: 
              fi.write(" %d %d %d\n" % (-1, nobs+1, -1) ) # First obs.
          elif nobs == data.size:
              fi.write(" %d %d %d\n" % (nobs-1, -1, -1) ) # Last obs.
          else:
              fi.write(" %d %d %d\n" % (nobs-1, -1, -1) ) 
      
          fi.write("obdef\n")
          fi.write("loc3d\n")

          fi.write("    %20.14f          %20.14f          %20.14f\n" % (lons[j,i], lats[j,i], hgts[k,j,i]))
      
          fi.write("     %d     \n" % vert_coord )
      
          fi.write("kind\n")
          
          if kind == ObType_LookUp("REFLECTIVITY") and data[k,j,i] <= 0.1 and zero_dbz_obtype:
              fi.write("     %d     \n" % ObType_LookUp("RADAR_CLEARAIR_REFLECTIVITY") )
              nobs_clearair += 1
          else:     
              fi.write("     %d     \n" % kind )

# If this GEOS cloud pressure observation, write out extra information (NOTE - NOT TESTED FOR HDF2ASCII LJW 04/13/15)
# 
#       if kind == ObType_LookUp("GOES_CWP_PATH"):
#           fi.write("    %20.14f          %20.14f  \n" % (row["satellite"][0], row["satellite"][1]) )
#           fi.write("    %20.14f  \n" % (row["satellite"][2]) )

# Check to see if its radial velocity and add platform informationp...need BETTER CHECK HERE!
      
          if kind == ObType_LookUp("VR"):
          
              R_xy            = np.sqrt(obs.xg[i]**2 + obs.yg[j]**2)
              elevation_angle = beam_elv(R_xy, obs.zg[k,j,i])

              platform_dir1 = (obs.xg[i] / R_xy) * np.deg2rad(elevation_angle)
              platform_dir2 = (obs.yg[j] / R_xy) * np.deg2rad(elevation_angle)
              platform_dir3 = np.sin(np.deg2rad(elevation_angle))
          
              fi.write("platform\n")
              fi.write("loc3d\n")

              fi.write("    %20.14f          %20.14f        %20.14f\n" % (platform_lon, platform_lat, platform_hgt) )
              fi.write("     %d     \n" % platform_vert_coord )
          
              fi.write("dir3d\n")
          
              fi.write("    %20.14f          %20.14f        %20.14f\n" % (platform_dir1, platform_dir2, platform_dir3) )
              fi.write("     %20.14f     \n" % obs.nyquist[k] )
              fi.write("     %d     \n" % platform_key )

    # Done with special radial velocity obs back to dumping out time, day, error variance info
      
          fi.write("    %d          %d     \n" % (seconds, days) )

    # Logic for command line override of observational error variances

          fi.write("    %20.14f  \n" % obs_error**2 )

          if nobs % 1000 == 0: print(" write_DART_ascii:  Processed observation # %d" % nobs)
  
      it.iternext()
      
  fi.close()
  
# To write out header information AFTER we know how big the observation data set is, we have
# to read back in the entire contents of the obs-seq file, store it, rewrite the file
# with header information first, and then dump the contents of obs-seq back inp.  Yuck.

  with file(filename, 'r') as f: f_obs_seq = f.read()

  fi = open(filename, "w")
  
  fi.write(" obs_sequence\n")
  fi.write("obs_kind_definitions\n")

# Deal with case that for reflectivity, 2 types of observations might have been created

  if kind == ObType_LookUp("REFLECTIVITY") and zero_dbz_obtype and nobs_clearair > 0:
      fi.write("       %d\n" % 2)
      akind, DART_name = ObType_LookUp(obs.field.upper(), DART_name=True)
      fi.write("    %d          %s   \n" % (akind, DART_name) )
      akind, DART_name = ObType_LookUp("RADAR_CLEARAIR_REFLECTIVITY", DART_name=True) 
      fi.write("    %d          %s   \n" % (akind, DART_name) )
  else:
      fi.write("       %d\n" % 1)
      akind, DART_name = ObType_LookUp(obs.field.upper(), DART_name=True)
      fi.write("    %d          %s   \n" % (akind, DART_name) )

  fi.write("  num_copies:            %d  num_qc:            %d\n" % (1, 1))
  
  fi.write(" num_obs:       %d  max_num_obs:       %d\n" % (nobs, nobs) )
      
  fi.write("observations\n")
  fi.write("QC radar\n")
          
  fi.write("  first:            %d  last:       %d\n" % (1, nobs) )

# Now write back in all the actual DART obs data

  fi.write(f_obs_seq)
  
  fi.close()
  
  print("\n write_DART_ascii:  Created ascii DART file, N = %d written" % nobs)
  
  if kind == ObType_LookUp("REFLECTIVITY") and zero_dbz_obtype and nobs_clearair > 0:
      print(" write_DART_ascii:  Number of clear air obs:             %d" % nobs_clearair)
      print(" write_DART_ascii:  Number of non-zero reflectivity obs: %d" % (nobs - nobs_clearair))

  return
  
########################################################################
# Volume Prep:  Threshold data and unfold velocities

def volume_prep(radar, unfold_type="phase"):
          
# Mask data beyond max_range

  max_range_gate = np.abs(radar.range['data'] - _radar_parameters['max_range']).argmin()
  radar.fields['reflectivity']['data'][:,max_range_gate:] = np.ma.masked
  radar.fields['velocity']['data'][:,max_range_gate:] = np.ma.masked

# Filter based on masking, dBZ threshold, and invalid gates

  gatefilter = pyart.correct.GateFilter(radar)
  gatefilter.exclude_invalid('velocity')
  gatefilter.exclude_invalid('reflectivity')
  gatefilter.exclude_masked('velocity')
  gatefilter.exclude_masked('reflectivity')
  gatefilter.exclude_below('reflectivity', _radar_parameters['min_dbz_analysis'])

# Dealias the velocity data

  if unfold_type == "phase":
    dealiased_radar = pyart.correct.dealias_unwrap_phase(radar, unwrap_unit='sweep', 
                                   nyquist_vel=None, check_nyquist_uniform=True, 
                                   gatefilter=gatefilter, rays_wrap_around=None, 
                                   keep_original=False, set_limits=True, 
                                   vel_field='velocity', corr_vel_field=None, 
                                   skip_checks=False)
                                   
    radar.add_field('unfolded velocity', dealiased_radar)
                                   

  if unfold_type == "region":
    dealiased_radar = pyart.correct.dealias_region_based(radar, interval_splits=3, 
                                   interval_limits=None, skip_between_rays=100, 
                                   skip_along_ray=100, centered=True, 
                                   nyquist_vel=None, check_nyquist_uniform=True, 
                                   gatefilter=False, rays_wrap_around=None, 
                                   keep_original=False, set_limits=True, 
                                   vel_field='velocity', corr_vel_field=None)
     
    radar.add_field('unfolded velocity', dealiased_radar)


#  Must implement function to get sounding data or specify previously unfolded radar 
#       volume to use PyART 4DD method. See PyART docs for more info.
#  if unfold_type == "4dd":  
#    snd_hgts,sngd_spds,snd_dirs = get_sound_data()
#    dealiased_radar = pyart.correct.dealias_fourdd(radar, last_radar=None, 
#                              sounding_heights=snd_hgts, sounding_wind_speeds=snd_spds, 
#                              sounding_wind_direction=snd_dirs, gatefilter=False, 
#                              filt=1, rsl_badval=131072.0, 
#                              keep_original=False, set_limits=True, 
#                              vel_field='velocity', corr_vel_field=None, 
#                              last_vel_field=None, debug=False, 
#                              max_shear=0.05, sign=1)

#   copies reflectivity data to a variable that doesn't exceed the 8 character variable name limitation in OPAWS 
#   radar.add_field_like('reflectivity', 'REF',
#                        radar.fields['reflectivity']['data'].copy())

  return gatefilter

########################################################################
#
# Grid data using parameters defined above in grid_dict 

def grid_data(volume, field):
  
  grid_spacing_xy = _grid_dict['grid_spacing_xy']
  grid_length     =_grid_dict['grid_radius_xy']
  grid_pts_xy     = 2*np.int(grid_length/grid_spacing_xy)
  roi             = _grid_dict['ROI']
  nthreads        = _grid_dict['nthreads']
  
  print '\n Gridding radar data with following parameters'
  print ' ---------------------------------------------\n'
  print ' Horizontal grid spacing: {} m'.format(grid_spacing_xy)
  print ' Grid points in x,y:      {},{}'.format(int(grid_pts_xy),int(grid_pts_xy))
  print ' Weighting function:      {}'.format(_grid_dict['weight_func'])
  print ' Radius of Influence:     {} km'.format(_grid_dict['ROI'])
  print ' Map projection:          {}'.format(_grid_dict['projection'])
  print ' Field to be gridded:     {}\n'.format(field) 
  print ' ---------------------------------------------\n' 
  
  area_id = 'Analysis grid'
  area_name = 'Analysis grid def'
  proj_id = 'lcc'
  radar_lat = volume.latitude['data'][0]
  radar_lon = volume.longitude['data'][0]
  proj4_args = '+proj=lcc +lat_0=%f +lon_0=%f +a=6371228.0 +units=m' % (radar_lat, radar_lon)  
  x_size = grid_pts_xy
  y_size = grid_pts_xy
  area_extent = (-grid_length, -grid_length, grid_length, grid_length)
  area_def    = utils.get_area_def(area_id, area_name, proj_id, proj4_args, x_size, y_size, area_extent)
  
  xg = area_def.proj_x_coords
  yg = area_def.proj_y_coords
  
  lons, lats = area_def.get_lonlats()
  
# Cressman weighting function

  wf = lambda r: (roi**2 - r**2) / (roi**2 + r**2)
  
# Create a 3D array for analysis grid, the vertical dimension is the number of tilts

  new = np.ma.zeros((volume.nsweeps, grid_pts_xy, grid_pts_xy))
  elevations = np.zeros((volume.nsweeps,))
  zgrid = np.zeros((volume.nsweeps, grid_pts_xy, grid_pts_xy))
  nyquist = np.zeros((volume.nsweeps,))

  tt = timeit.clock()
  
  for n, refl_sweep_data in enumerate(volume.iter_field(field)):
    if n == 0:
      begin = 0
      end   = volume.sweep_end_ray_index['data'][n] + 1
    else:
      begin = volume.sweep_end_ray_index['data'][n-1] + 1
      end   = volume.sweep_end_ray_index['data'][n] + 1

    yob           = volume.gate_latitude['data'][begin:end,:]
    xob           = volume.gate_longitude['data'][begin:end,:]
    zob           = volume.gate_z['data'][begin:end,:]
    elevations[n] = volume.elevation['data'][begin:end].mean()
    nyquist[n]    = volume.get_nyquist_vel(n)
    
    x, y, z = volume.get_gate_x_y_z(n)
    
    obs_def = geometry.SwathDefinition(lons=xob, lats=yob) 

# Original method - but the object tree (obs_def) is created twice.  So this is a bit faster
#
#     tmp = kd_tree.resample_custom(obs_def, refl_sweep_data, \
#                                        area_def, radius_of_influence=roi, reduce_data=True, \
#                                        weight_funcs=wf, fill_value=np.nan, nprocs=nthreads)

    input_index, output_index, index_array, distances = \
          kd_tree.get_neighbour_info(obs_def, area_def, roi)
                                       
    tmp = kd_tree.get_sample_from_neighbour_info('custom', area_def.shape, refl_sweep_data,  \
                                                  input_index, output_index, index_array, \
                                                  distance_array=distances, weight_funcs=wf, fill_value=np.nan)

    if field == "reflectivity":
      new[n,...] = np.ma.masked_where(tmp < _radar_parameters['min_dbz_analysis'], tmp)
    else:
      new[n,...] = tmp
      
    zgrid[n,...] = kd_tree.get_sample_from_neighbour_info('custom', area_def.shape, z,  \
                                                  input_index, output_index, index_array, \
                                                  distance_array=distances, weight_funcs=wf, fill_value=np.nan)

  print("\n %f secs to run pyresample analysis for all levels \n" % (timeit.clock()-tt))

  return Gridded_Field("data_grid", field = field, data = new, proj4 = proj4_args, 
                       xg = xg, yg = yg, zg = zgrid,                   
                       lats = lats, lons = lons, elevations=elevations,
                       radar_lat = radar_lat, radar_lon=radar_lon, radar_hgt=volume.altitude['data'][0],
                       time = volume.time, metadata = volume.metadata, nyquist = nyquist  ) 

###########################################################################################
#
# Read environment variables for plotting a shapefile

def plot_shapefiles(map, shapefiles=None, color='k', linewidth=0.5, counties=False, ax=None):

    if shapefiles:
        try:
            shapelist = os.getenv(shapefiles).split(":")

            if len(shapelist) > 0:

                for item in shapelist:
                    items      = item.split(",")
                    shapefile  = items[0]
                    color      = items[1]
                    linewidth  = items[2]

                    s = map.readshapefile(shapefile,'myshapes',drawbounds=False)

                    for shape in map.counties:
                        xx, yy = zip(*shape)
                        map.plot(xx, yy, color=color, linewidth=linewidth, ax=ax, zorder=4)

        except OSError:
            print "PLOT_SHAPEFILES:  NO SHAPEFILE ENV VARIABLE FOUND "
            
    if counties:
            map.drawcounties(ax=ax, linewidth=0.5, color='k', zorder=5)
            
########################################################################
#
# Create two panel plot of processed, gridded velocity and reflectivity data  

def grid_plot(ref, vel, sweep, fsuffix=None, shapefiles=None, interactive=True):
  
# Set up colormaps 

  from matplotlib.colors import BoundaryNorm
   
  cmapr = cm.NWSRef
  cmapr.set_bad('white',1.0)
  cmapr.set_under('white',1.0)

  cmapv = cm.BuDRd18
  cmapv.set_bad('white',1.)
  cmapv.set_under('black',1.)
  cmapv.set_over('black',1.)
  
  normr = BoundaryNorm(np.arange(10, 85, 5), cmapr.N)
  normv = BoundaryNorm(np.arange(-32, 34, 2), cmapv.N)
  
  min_dbz = _radar_parameters['min_dbz_analysis']  
  xwidth = ref.xg.max() - ref.xg.min()
  ywidth = ref.yg.max() - ref.yg.min()

# Create png file label

  if fsuffix == None:
      print("\n pyROTH.grid_plot:  No output file name is given, writing to %s" % "VR_RF_...png")
      filename = "VR_RF_%2.2d_plot.png" % (sweep)
  else:
       filename = "VR_RF_%2.2d_%s.png" % (sweep, fsuffix)

  fig, (ax1, ax2) = P.subplots(1, 2, sharey=True, figsize=(15,8))
  
# Set up coordinates for the plots
  bgmap = Basemap(projection=_grid_dict['projection'], width=xwidth, \
                  height=ywidth, resolution='c', lat_0=ref.radar_lat,lon_0=ref.radar_lon, ax=ax1)
                  
  xoffset,yoffset = bgmap(ref.radar_lon,ref.radar_lat)
  xg = xoffset + ref.xg
  yg = yoffset + ref.yg
  
  yg_2d, xg_2d = np.meshgrid(ref.yg, ref.xg)
  
  yg_2d = yg_2d.transpose() + yoffset
  xg_2d = xg_2d.transpose() + xoffset
  
# fix xg, yg coordinates so that pcolormesh plots them in the center.

  dx2 = 0.5*(ref.xg[1] - ref.xg[0])
  dy2 = 0.5*(ref.yg[1] - ref.yg[0])
  
  xe = np.append(xg-dx2, [xg[-1] + dx2])
  ye = np.append(yg-dy2, [yg[-1] + dy2])

# REFLECTVITY PLOT

  if shapefiles:
      plot_shapefiles(bgmap, shapefiles=shapefiles, counties=_plot_counties, ax=ax1)
  else:
      plot_shapefiles(bgmap, counties=_plot_counties, ax=ax1)
 
  bgmap.drawparallels(range(10,80,1),    labels=[1,0,0,0], linewidth=0.5, ax=ax1)
  bgmap.drawmeridians(range(-170,-10,1), labels=[0,0,0,1], linewidth=0.5, ax=ax1)

  im1 = bgmap.pcolormesh(xe, ye, ref.data[sweep], cmap=cmapr, norm=normr, ax=ax1)
  cbar = bgmap.colorbar(im1,location='right')
  cbar.set_label('Reflectivity (dBZ)')
  ax1.set_title('Thresholded Reflectivity (Gridded)')
  bgmap.scatter(xoffset,yoffset, c='k', s=50., alpha=0.8, ax=ax1)

# Plot missing values as points.
#   r_mask = (ref.data.mask[sweep] == True)
#   bgmap.scatter(xg_2d[r_mask].transpose(), yg_2d[r_mask].transpose(), c='k', s = 1., alpha=0.5, ax=ax1)

# Plot zeros as "o"
  r_mask = np.logical_and(ref.data[sweep] < 1.0, (ref.data.mask[sweep] == False))
  bgmap.scatter(xg_2d[r_mask].transpose(), yg_2d[r_mask].transpose(), s=25, facecolors='none', edgecolors='k', alpha=1.0, ax=ax1)
  
# RADIAL VELOCITY PLOT

  bgmap = Basemap(projection=_grid_dict['projection'], width=xwidth, \
                  height=ywidth, resolution='c', lat_0=ref.radar_lat,lon_0=ref.radar_lon, ax=ax2)
                  
  if shapefiles:
      plot_shapefiles(bgmap, shapefiles=shapefiles, counties=_plot_counties, ax=ax2)
  else:
      plot_shapefiles(bgmap, counties=_plot_counties, ax=ax2)
    
  bgmap.drawparallels(range(10,80,1),labels=[1,0,0,0], linewidth=0.5, ax=ax2)
  bgmap.drawmeridians(range(-170,-10,1),labels=[0,0,0,1],linewidth=0.5, ax=ax2)
  
  im1 = bgmap.pcolormesh(xe, ye, vel.data[sweep], cmap=cmapv, norm=normv, ax=ax2)
  cbar = bgmap.colorbar(im1,location='right')
  cbar.set_label('Dealised Radial Velocity (meters_per_second)')
  ax2.set_title('Thresholded, Unfolded Radial Velocity (Gridded)') 
  bgmap.scatter(xoffset,yoffset, c='k', s=50.,alpha=0.8, ax=ax2)
  
# Now plot locations of nan data

  v_mask = vel.data.mask == True
  bgmap.scatter(xg_2d[v_mask[sweep]], yg_2d[v_mask[sweep]], c='k', s = 1., alpha=0.5, ax=ax2)

# Get other metadata....for labeling
  instrument_name = ref.metadata['instrument_name']
  time_start = ncdf.num2date(ref.time['data'][0], ref.time['units'])
  time_text = time_start.isoformat().replace("T"," ")
  title = '\nDate:  %s   Time:  %s Z   Elevation:  %2.2f deg' % (time_text[0:10], time_text[10:19], ref.elevations[sweep])
  P.suptitle(title, fontsize=24)
  
  P.savefig(filename)
  
  if interactive:  P.show()
  
########################################################################
# Main function

if __name__ == "__main__":

  print ' ================================================================================'
  print ''
  print '                   BEGIN PROGRAM pyROTH                     '
  print ''

  parser = OptionParser()
  parser.add_option("-f", "--file",       dest="fname",     default=None,  type="string", \
                     help = "filename of NEXRAD level II volume to process")
  parser.add_option("-u", "--unfold",     dest="unfold",    default="phase",  type="string", \
                     help = "dealiasing method to use (phase or region, default = phase)")
  parser.add_option("-w", "--write",     dest="write",   default=False, \
                     help = "Boolean flag to write DART ascii file", action="store_true")
  parser.add_option("-p", "--plot",       dest="plot",      default=-1,  type="int",      \
                     help = "Specify a number between 0 and # elevations to plot ref and vr in that co-plane")
  parser.add_option("-i", "--interactive", dest="interactive", default=False,  action="store_true",     \
                     help = "Boolean flag to specify to plot image to screen (when plot > -1).")  
  parser.add_option("-s", "--shapefiles", dest="shapefiles", default=None, type="string",    \
                     help = "Name of system env shapefile you want to add to the plots.")
                     
                   
  (options, args) = parser.parse_args()
  
  parser.print_help()

  print ''
  print ' ================================================================================'

  if options.fname == None:
    print "\n\n ***** USER MUST SPECIFY NEXRAD LEVEL II (MESSAGE 31) FILE! *****"
    print "\n                         EXITING!\n\n"
    parser.print_help()
    print
    sys.exit(1)
  else:
    fname   = os.path.abspath(options.fname)
    fsuffix = os.path.split(fname)[-1][0:17]
    fsuffix = fsuffix[0:4] + "_" + fsuffix[4:]
  
  if options.unfold == "phase":
    print "\n pyROTH dealias_unwrap_phase unfolding will be used\n"
    unfold_type = "phase"
  elif options.unfold == "region":
    print "\n pyROTH dealias_region_based unfolding will be used\n"
    unfold_type = "region"
  else:
    print "\n ***** INVALID OR NO VELOCITY DEALIASING METHOD SPECIFIED *****"
    print "\n          NO VELOCITY UNFOLDING DONE...\n\n"
    unfold_type = None

  if options.plot < 0:
    plot_grid = False
  else:
    sweep_num = options.plot
    plot_grid = True
  
# Read input file and create radar object

  print '\n Reading: {}\n'.format(fname)
   
  t0   = timeit.time()
   
  tim0 = timeit.time()
  
  volume = pyart.io.read_nexrad_archive(fname, field_names=None, 
                                        additional_metadata=None, file_field_names=False, 
                                        delay_field_loading=False, 
                                        station=None, scans=None, linear_interp=True)
  pyROTH_io_cpu = timeit.time() - tim0
  
  print "\n Time for reading in LVL2: {} seconds".format(pyROTH_io_cpu)
  
  tim0 = timeit.time()
  
  gatefilter = volume_prep(volume, unfold_type=unfold_type) 
  
  pyROTH_unfold_cpu = timeit.time() - tim0

  print "\n Time for unfolding velocity: {} seconds".format(pyROTH_unfold_cpu)
  
  tim0 = timeit.time()

  ref = dbz_masking(grid_data(volume, "reflectivity"))

  if unfold_type == None:  
      vel = grid_data(volume, "velocity")
  else:
      vel = grid_data(volume, "unfolded velocity")

  pyROTH_regrid_cpu = timeit.time() - tim0
  
  print "\n Time for gridding fields: {} seconds".format(pyROTH_regrid_cpu)
    
  if plot_grid:
      plottime = grid_plot(ref, vel, sweep_num, fsuffix=fsuffix, shapefiles=options.shapefiles, interactive=options.interactive)

  if options.write == True:      
      ret = write_DART_ascii(vel, fsuffix=fsuffix+"_VR", obs_error=_obs_errors['velocity'])
      ret = write_DART_ascii(ref, fsuffix=fsuffix+"_RF", obs_error=_obs_errors['reflectivity'])
    
  pyROTH_cpu_time = timeit.time() - t0
  
  print "\n Time for pyROTH operations: {} seconds".format(pyROTH_cpu_time)

  print "\n PROGRAM pyROTH COMPLETED\n"
