#!/usr/bin/env python
#############################################################
# grid3d: A program to process new MRMS volumes             #
#         netCDF4 MRMS files are read in and processed to   #
#         to DART format for assimilationp.                 #
#                                                           #
#       Python package requirements:                        #
#       ----------------------------                        #
#       Numpy                                               #
#       Scipy                                               #
#       matplotlib                                          #
#############################################################
#
# created by Lou Wicker Feb 2017
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
from optparse import OptionParser
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText
import netCDF4 as ncdf

from pyproj import Proj
import pylab as plt  
from mpl_toolkits.basemap import Basemap
from pyart.graph import cm
import datetime as DT
from numpy import ma
from dart_tools import *

# missing value
_missing = -9999.
_plot_counties = True

# Colorscale information
_ref_scale = (0.,74.)
_vr_scale  = (-40.,40.)

# Colortables
_ref_ctable = cm.NWSRef
_vr_ctable  = cm.Carbone42

##########################################################################################
# Parameter dict for reflectivity masking
_grid_dict = {
              'zero_dbz_obtype' : True,
              'thin_zeros'      : 3,
              'halo_footprint'  : 4,
              'max_height'      : 10000.,
              'zero_levels'     : [5000.], 
              'min_dbz_analysis': 15.0,
              'min_dbz_zeros'   : 15.0,
              'reflectivity'    : 5.0,
              '0reflectivity'   : 5.0, 
              'levels'          : [1,2,3,4,5,6,8,10,12,14, 16, 18],
              'QC_info'         : [[15.,5.],[20.,1.]],
             }

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
# Get the filenames out of the directory

def get_dir_files(dir, pattern, Quiet=False):

    try:
        os.path.isdir(dir)
    except:
        print(" %s is not a directory, exiting" % dir)
        sys.exit(1)
        
    filenames = glob.glob("%s/Reflect*.nc" % (dir))

    if not Quiet:
        print("\n Processing %d files in the directory:  %s" % (len(filenames), dir))
        print("\n First file is %s" % (os.path.basename(filenames[0])))
        print("\n Last  file is %s\n" % (os.path.basename(filenames[-1])))

    return filenames

#=========================================================================================
# DBZ Mask

def dbz_masking(ref, thin_zeros=2):

   nz, ny, nx = ref.data.shape
  
# First create zeros
 
   zeros = np.ones((ny, nx), dtype=np.float32)

   zero_dbz = np.ma.masked_where(zeros == 1.0, zeros)
   
   print(" Initial number of zeros in the field:  %d\n " % (np.sum(zero_dbz.mask==False)))
   
   if thin_zeros > 0:
      zero_dbz.data[::thin_zeros, ::thin_zeros] = 0.0       # Here is where the thinned number of zeros is created
      zero_dbz.mask[::thin_zeros, ::thin_zeros] = False     # because where mask == False, we have a zero...

   print(" Number of zeros in the field after thining:  %d\n " % (np.sum(zero_dbz.mask==False)))
   
   ref.composite = ref.data.data.max(axis=0)
   
   zero_dbz.mask[ ref.composite > _grid_dict['min_dbz_zeros'] ] = True

   print(" Number of zeros in the field after compositive ref mask:  %d\n " % (np.sum(zero_dbz.mask==False)))
      
   raw_field = np.where(ref.composite < _grid_dict['min_dbz_zeros'], 0.0, ref.composite)
     
   max_neighbor = (ndimage.maximum_filter(raw_field, size=_grid_dict['halo_footprint']) > 0.1)
     
   zero_dbz.mask[max_neighbor == True] = True

   print(" Number of zeros in the field after halo:  %d\n " % (np.sum(zero_dbz.mask==False)))
   
# Create new object so that we can combine later

   obj_zero_dbz = Gridded_Field("zero_dbz", data=zero_dbz, zg=_grid_dict['zero_levels'], radar_hgt=0.0)
   
   ref.zero_dbz = obj_zero_dbz
   
# Now clip the reflectivity

   ref.data.mask[ ref.data < _grid_dict['min_dbz_analysis'] ] = True
             
   if _grid_dict['max_height'] > 0:
       zg3d          = np.repeat(ref.zg, ny*nx).reshape(nz,ny,nx)
       mask1         = (zg3d > _grid_dict['max_height'])  
       mask2         = ref.data.mask
       ref.data.mask = np.logical_or(mask1, mask2)
        
   return ref

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

def grid_plot(ref, sweep, fsuffix=None, shapefiles=None, interactive=True):
  
# Set up colormaps 

  from matplotlib.colors import BoundaryNorm
   
  cmapr = cm.NWSRef
  cmapr.set_bad('white',1.0)
  cmapr.set_under('white',1.0)

  cmapv = cm.Carbone42
  cmapv.set_bad('white',1.)
  cmapv.set_under('black',1.)
  cmapv.set_over('black',1.)
  
  normr = BoundaryNorm(np.arange(10, 85, 5), cmapr.N)
  normv = BoundaryNorm(np.arange(-48, 50, 2), cmapv.N)
  
# Create png file label

  if fsuffix == None:
      print("\n pyROTH.grid_plot:  No output file name is given, writing to %s" % "RF_...png")
      filename = "RF_%2.2d_plot.png" % (sweep)
  else:
       filename = "RF_%2.2d_%s.png" % (sweep, fsuffix.split("/")[1])

  fig, ax1 = plt.subplots(figsize=(14,10))
  
# Set up coordinates for the plots

  sw_lon = ref.lons.min()
  ne_lon = ref.lons.max()
  sw_lat = ref.lats.min()
  ne_lat = ref.lats.max()

  bgmap = Basemap(projection='lcc', llcrnrlon=sw_lon,llcrnrlat=sw_lat,urcrnrlon=ne_lon,urcrnrlat=ne_lat, \
              lat_0=0.5*(ne_lat+sw_lat), lon_0=0.5*(ne_lon+sw_lon), resolution='c', ax=ax1)
                  
  xg, yg = bgmap(ref.lons, ref.lats)
    
  yg_2d, xg_2d = np.meshgrid(yg, xg)
  
  yg_2d = yg_2d.transpose()
  xg_2d = xg_2d.transpose()
  
# fix xg, yg coordinates so that pcolormesh plots them in the center.

  dx2 = 0.5*(xg[1] - xg[0])
  dy2 = 0.5*(yg[1] - yg[0])
  
  xe = np.append(xg-dx2, [xg[-1] + dx2])
  ye = np.append(yg-dy2, [yg[-1] + dy2])


# REFLECTVITY PLOT

  if shapefiles:
      plot_shapefiles(bgmap, shapefiles=shapefiles, counties=_plot_counties, ax=ax1)
  else:
      plot_shapefiles(bgmap, counties=_plot_counties, ax=ax1)
 
  bgmap.drawparallels(range(10,80,1),    labels=[1,0,0,0], linewidth=0.5, ax=ax1)
  bgmap.drawmeridians(range(-170,-10,1), labels=[0,0,0,1], linewidth=0.5, ax=ax1)

  im1 = bgmap.pcolormesh(xe, ye, ref.data[sweep], cmap=cmapr, vmin = _ref_scale[0], vmax = _ref_scale[1], ax=ax1)
  cbar = bgmap.colorbar(im1,location='right')
  cbar.set_label('Reflectivity (dBZ)')
  ax1.set_title('Thresholded Reflectivity (Gridded)')
  
  at = AnchoredText("Max dBZ: %4.1f" % (ref.data[sweep].max()), loc=4, prop=dict(size=12), frameon=True,)
  at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
  ax1.add_artist(at)

# Plot missing values as points.
#   r_mask = (ref.data.mask[sweep] == True)
#   bgmap.scatter(xg_2d[r_mask], yg_2d[r_mask], c='k', s = 1., alpha=0.5, ax=ax1)

# Plot zeros as "o"

  try:
      r_mask = (ref.zero_dbz.mask == False)
      print("%d" % np.sum(r_mask) )
      bgmap.scatter(xg_2d[r_mask], yg_2d[r_mask], s=25, facecolors='none', \
                    edgecolors='k', alpha=1.0, ax=ax1) 
  except AttributeError:
      r_mask = (ref.data.data[sweep] < 1.0) & (ref.data.mask[sweep] == False)
      print("%d" % np.sum(r_mask) )
      bgmap.scatter(xg_2d[r_mask], yg_2d[r_mask], s=25, facecolors='none', \
                    edgecolors='k', alpha=1.0, ax=ax1)
  
# Get other metadata....for labeling

#   time_start = ncdf.num2date(ref.time['data'][0], ref.time['units'])
  time_text = ref.time.strftime('%Y-%m-%d %H:%M')
  title = '\nDate:  %s   Time:  %s Z   Height:  %4.1f km' % (time_text[0:10], time_text[10:19], 0.001*ref.zg[sweep])
  plt.suptitle(title, fontsize=24)
  
  plt.savefig(filename)
  
  if interactive:  plt.show()
#-------------------------------------------------------------------------------
# Main function defined to return correct sys.exit() calls

def main(argv=None):
   if argv is None:
       argv = sys.argv
#
# Command line interface 
#
   parser = OptionParser()
   parser.add_option("-d", "--dir",   dest="dir",    default=None,  type="string", help = "Directory where MRMS files are")
   
   parser.add_option(      "--realtime",  dest="realtime",   default=None,  \
               help = "Boolean flag to uses this YYYYMMDDHHMM time stamp for the realtime processing")
 
   parser.add_option("-g", "--grep",  dest="grep",   default="*.nc", type="string", help = "Pattern grep, [*.nc, *VR.h5]")
   
   parser.add_option("-w", "--write", dest="write",   default=False, \
                           help = "Boolean flag to write DART ascii file", action="store_true")
                           
   parser.add_option("-o", "--out",      dest="out_dir",  default="ref_files",  type="string", \
                           help = "Directory to place output files in")
                           
   parser.add_option("-p", "--plot",      dest="plot",      default=-1,  type="int",      \
                     help = "Specify a number between 0 and 20 to plot reflectivity")
                  
   (options, args) = parser.parse_args()

#-------------------------------------------------------------------------------

   if options.dir == None:
          
      print "\n\n ***** USER MUST SPECIFY A DIRECTORY WHERE FILES ARE *****"
      print "\n                         EXITING!\n\n"
      parser.print_help()
      print
      sys.exit(1)
      
   if options.plot < 0:
       plot_grid = False
   else:
       sweep_num = options.plot
       plot_grid = True

   if options.realtime != None:
       year = int(options.realtime[0:4])
       mon  = int(options.realtime[4:6])
       day  = int(options.realtime[6:8])
       hour = int(options.realtime[8:10])
       minute = int(options.realtime[10:12])
       a_time = DT.datetime(year, mon, day, hour, minute, 0)

#-------------------------------------------------------------------------------

# if realtime, now find the file created within T+5 of analysis

   if options.realtime != None:

       in_filenames = get_dir_files(options.dir, options.grep, Quiet=True)

       for n, file in enumerate(in_filenames):
           f_time = DT.datetime.strptime(os.path.basename(file)[-18:-3], "%Y%m%d-%H%M%S")
           if (f_time > a_time - DT.timedelta(minutes = 1)) and (f_time <= a_time + DT.timedelta(minutes = 6)):
                last_file    = in_filenames[n]

       in_filenames = [last_file]
       print("\n prep_grid3d:  RealTime FLAG is true, only processing %s\n" % (in_filenames[:]))
       rlt_filename = "%s_%s" % ("obs_seq_RF", a_time.strftime("%Y%m%d%H%M"))

   else:
       in_filenames = get_dir_files(options.dir, options.grep, Quiet=False)
       out_filenames = []
       print("\n prep_grid3d:  Processing %d files in the directory:  %s\n" % (len(in_filenames), options.dir))
       print("\n prep_grid3d:  First file is %s\n" % (in_filenames[0]))
       print("\n prep_grid3d:  Last  file is %s\n" % (in_filenames[-1]))
   
 # Make sure there is a directory to write files into....
 
   if not os.path.exists(options.out_dir):
       try:
           os.mkdir(options.out_dir)
       except:
           print("\n**********************   FATAL ERROR!!  ************************************")
           print("\n PREP_GRID3D:  Cannot create output dir:  %s\n" % options.out_dir)
           print("\n**********************   FATAL ERROR!!  ************************************")

   for n, file in enumerate(in_filenames):
   
      print ' ================================================================================'

      if options.realtime != None:
          out_filename = os.path.join(options.out_dir, rlt_filename)
          time         = a_time
	  print(" Out filename:  %s\n" % out_filename)
      else:
          print options.realtime == None
          str_time     = "%s_%s" % (os.path.basename(file)[-18:-10], os.path.basename(file)[-9:-3])
          prefix       = "obs_seq_RF_%s" % str_time
          out_filename = os.path.join(options.out_dir, prefix)
          time         = DT.datetime.strptime(file[-18:-3], "%Y%m%d-%H%M%S")
          print(" Out filename:  %s\n" % out_filename)
   
      f       = ncdf.Dataset(file, "r")
      missing = f.MissingData
      nlon    = len(f.dimensions['Lon'])
      nlat    = len(f.dimensions['Lat'])
      nlvl    = len(f.dimensions['Ht'])
      ref     = ma.MaskedArray(f.variables['ReflectivityQC'][...], mask = (f.variables['ReflectivityQC'][...] < missing+1.))
      lons    = f.variables['Lon'][...]
      lats    = f.variables['Lat'][...]
      msl     = f.variables['Height'][...]
      height0 = f.Height
      dtime   = (time - DT.datetime(1970, 1, 1, 0, 0)).total_seconds

      f.close()
   
      print(" Time of file:  %s\n" % time.strftime('%Y-%m-%d %H:%M:%S') )
             
      ref_obj = Gridded_Field(file, data = ref, field = "REFLECTIVITY", zg = msl, \
                                    lats = lats, lons = lons, radar_hgt = height0, time = time ) 
#       ref_obj.data.field     = "REFLECTIVITY"
#       ref_obj.data.zg        = msl
#       ref_obj.data.lats      = lats
#       ref_obj.data.lons      = lons
#       ref_obj.data.radar_hgt = height0
# #       ref_obj.data.time      = {'data': np.array(dtime), 'units': "seconds since January 1, 1970"}
#       ref_obj.data.time      = time
      
      ref_obj = dbz_masking(ref_obj, thin_zeros=_grid_dict['thin_zeros'])
      
      if plot_grid:
          grid_plot(ref_obj, sweep_num)
    
      if options.write == True:      
          ret = write_DART_ascii(ref_obj, filename=out_filename, levels=_grid_dict['levels'],
                                 obs_error=[_grid_dict['reflectivity'],_grid_dict['0reflectivity']], 
                                 QC_info=_grid_dict['QC_info'])


#-------------------------------------------------------------------------------
# Main program for testing...
#
if __name__ == "__main__":
    sys.exit(main())
