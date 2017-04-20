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

import warnings
warnings.filterwarnings("ignore")

# missing value
_missing = -9999.
_plot_counties = True

# Colorscale information

_ref_scale  = np.arange(5.,85.,5.0)
_ref_ctable = cm.NWSRef
_ref_min_plot = 20.
_plot_format = 'png'

# Radar information

_dbz_name         = "MergedReflectivityQC"
_temp_netcdf_file = "/work/wicker/REALTIME/tmp.nc"

# Grid stuff

NX = 180
NY = 180

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
              'levels'          : ['00.50','01.00','01.50','02.00','02.50','03.00','03.50', \
                                   '04.00','05.00','06.00','07.00','08.00','09.00','10.00'],
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

#===============================================================================
def get_loc(x, xc, radius):

  """
      Finds the array index of locations around a center point
  """
  if xc < 0:
      indices = np.where(x >= xc)
  else:
      indices = np.where(x < xc)

  if np.size(indices[0]) == 0:
    return -1, 0
  else:
    i0 = indices[0][0]
    return i0, i0+1
                 
#=========================================================================================
# Get the filenames out of the directory

def get_dir_files(dir, pattern, Quiet=False):

    try:
        os.path.isdir(dir)
    except:
        print(" %s is not a directory, exiting" % dir)
        sys.exit(1)
        
    filenames = glob.glob("%s/01.00/*.netcdf.gz" % (dir))

    if not Quiet:
        print("\n Processing %d files in the directory:  %s" % (len(filenames), dir))
        print("\n First file is %s" % (os.path.basename(filenames[0])))
        print("\n Last  file is %s\n" % (os.path.basename(filenames[-1])))

    return filenames

#=========================================================================================
# Get the filenames out of the directory

def assemble_3D_grid(path, filename, loc=None, Quiet=False):

    try:
        os.path.isdir(path)
    except:
        print(" %s is not a directory, exiting" % path)
        sys.exit(1)
        
    levels = _grid_dict['levels']
    nlvls  = len(levels)
    
    for n, lvl in enumerate(levels):
    
         full_path = os.path.join(path, lvl, filename)
         print full_path
         
         if not Quiet:
             print("\n Processing level %s in the directory:  %s" % (lvl, path))
      
         cmd = "gunzip -c %s > %s" % (full_path, _temp_netcdf_file)
         if not Quiet:
             print("\n Running......%s" % (cmd))

         os.system(cmd)
         
         f = ncdf.Dataset(_temp_netcdf_file, "r")
         
         if n == 0:
             
             nlons  = len(f.dimensions['Lon'])
             nlats  = len(f.dimensions['Lat'])
             f_lats = f.variables['Lat'][...]
             f_lons = f.variables['Lon'][...]
             
             missingData = f.MissingData
             time   = DT.datetime.fromtimestamp(f.variables['time'][0])
             
             if loc != None:
                 ic     = get_loc(f_lons, loc[1], 0.5)[0]
                 jc     = get_loc(f_lats, loc[0], 0.5)[0]
             else:
                 ic, jc = nlats/2, nlons/2

             i0, i1 = ic-NX/2, ic+NX/2
             j0, j1 = jc-NY/2, jc+NY/2

             if not Quiet:
                 print("\n Lon indexes:  %d  %d" % (i0, i1))
                 print("\n Lat indexes:  %d  %d" % (j0, j1))
                      
             g_lats    = f_lats[j0:j1]
             g_lons    = f_lons[i0:i1]
             array     = missingData * np.ones((nlvls, g_lats.size, g_lons.size))
             g_heights = np.zeros((nlvls,))

         g_heights[n] = f.Height    
         array[n,...] = f.variables[_dbz_name][0,i0:i1,j0:j1]
         
         f.close()
  
    ref = ma.MaskedArray(array, mask = (array < missingData+1.))        
    
    return Gridded_Field(filename, data = ref, field = "REFLECTIVITY", zg = g_heights, \
                         lats = g_lats, lons = g_lons, radar_hgt = 0.0, time = time, \
                         missingData = missingData ) 

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

def plot_shapefiles(map, shapefiles=None, color='k', linewidth=0.5, states=True, counties=False, ax=None):

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
            map.drawcounties(ax=ax, linewidth=0.25, color='0.5', zorder=5)

    if states:
            map.drawstates(ax=ax, linewidth=0.25, color='k', zorder=5)
            
########################################################################
#
# Create processed reflectivity data  

def grid_plot(ref, sweep, plot_filename=None, shapefiles=None, interactive=False):
  
# Set up colormaps 

  from matplotlib.colors import BoundaryNorm
   
  cmapr = _ref_ctable
  cmapr.set_bad('white',1.0)
  cmapr.set_under('white',1.0)

  normr = BoundaryNorm(_ref_scale, cmapr.N)
  
# Create png file label

  if plot_filename == None:
      print("\n pyROTH.grid_plot:  No output file name is given, writing to %s" % "RF_...png")
      filename = "RF_%2.2d_plot.%s" % (sweep, _plot_format)
  else:
       filename = "%s.%s" % (plot_filename, _plot_format)

  fig, axes = plt.subplots(1, 2, sharey=True, figsize=(15,10))
  
# Set up coordinates for the plots

  sw_lon = ref.lons.min()
  ne_lon = ref.lons.max()
  sw_lat = ref.lats.min()
  ne_lat = ref.lats.max()

  bgmap = Basemap(projection='lcc', llcrnrlon=sw_lon,llcrnrlat=sw_lat,urcrnrlon=ne_lon,urcrnrlat=ne_lat, \
                  lat_0=0.5*(ne_lat+sw_lat), lon_0=0.5*(ne_lon+sw_lon), resolution='c', ax=axes[0])

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
      plot_shapefiles(bgmap, shapefiles=shapefiles, counties=_plot_counties, ax=axes[0])
  else:
      plot_shapefiles(bgmap, counties=_plot_counties, ax=axes[0])
 
  bgmap.drawparallels(range(10,80,1),    labels=[1,0,0,0], linewidth=0.5, ax=axes[0])
  bgmap.drawmeridians(range(-170,-10,1), labels=[0,0,0,1], linewidth=0.5, ax=axes[0])

  if sweep == -1:
      ref_data = ref.data.max(axis=0)
  else:
      ref_data = ref.data[sweep]

# pixelated plot

  im1  = bgmap.pcolormesh(xe, ye, ref_data, cmap=cmapr, norm=normr, vmin = _ref_min_plot, vmax = _ref_scale.max(), ax=axes[0])
  cbar = bgmap.colorbar(im1,location='right')
  cbar.set_label('Reflectivity (dBZ)')
  axes[0].set_title('Pixel Reflectivity (Gridded)')
  at = AnchoredText("Max dBZ: %4.1f" % (ref_data.max()), loc=4, prop=dict(size=10), frameon=True,)
  at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
  axes[0].add_artist(at)

# Contour filled plot

  bgmap = Basemap(projection='lcc', llcrnrlon=sw_lon,llcrnrlat=sw_lat,urcrnrlon=ne_lon,urcrnrlat=ne_lat, \
                  lat_0=0.5*(ne_lat+sw_lat), lon_0=0.5*(ne_lon+sw_lon), resolution='c', ax=axes[1])

  if shapefiles:
      plot_shapefiles(bgmap, shapefiles=shapefiles, counties=_plot_counties, ax=axes[1])
  else:
      plot_shapefiles(bgmap, counties=_plot_counties, ax=axes[1])
 
  bgmap.drawparallels(range(10,80,1),    labels=[1,0,0,0], linewidth=0.5, ax=axes[1])
  bgmap.drawmeridians(range(-170,-10,1), labels=[0,0,0,1], linewidth=0.5, ax=axes[1])
                  
  im1 = bgmap.contourf(xg_2d, yg_2d, ref_data, levels= _ref_scale, cmap=cmapr, norm=normr, \
                       vmin = _ref_min_plot, vmax = _ref_scale.max(), ax=axes[1])

  cbar = bgmap.colorbar(im1,location='right')
  cbar.set_label('Reflectivity (dBZ)')
  axes[1].set_title('Contoured Reflectivity (Gridded)')
  at = AnchoredText("Max dBZ: %4.1f" % (ref_data.max()), loc=4, prop=dict(size=10), frameon=True,)
  at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
  axes[1].add_artist(at)

# Plot missing values as points.
#   r_mask = (ref.data.mask[sweep] == True)
#   bgmap.scatter(xg_2d[r_mask], yg_2d[r_mask], c='k', s = 1., alpha=0.5, ax=ax1)

# Plot zeros as "o"

  try:
      r_mask = (ref.zero_dbz.mask == False)
      print("%d" % np.sum(r_mask) )
      bgmap.scatter(xg_2d[r_mask], yg_2d[r_mask], s=25, facecolors='none', \
                    edgecolors='k', alpha=1.0, ax=axes[0]) 
  except AttributeError:
      print "Attribute error"
      r_mask = (ref.data.data[4] < 1.0) & (ref.data.mask[4] == False)
      print("%d" % np.sum(r_mask) )
      bgmap.scatter(xg_2d[r_mask], yg_2d[r_mask], s=25, facecolors='none', \
                    edgecolors='k', alpha=1.0, ax=axes[0])

# Get other metadata....for labeling

  time_text = ref.time.strftime('%Y-%m-%d %H:%M')

  title = '\nDate:  %s   Time:  %s Z   Height:  %4.1f km' % (time_text[0:10], time_text[10:19], 0.001*ref.zg[0])
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
 
   parser.add_option("-g", "--grep",  dest="grep",   default="*.netcdf.gz", type="string", help = "Pattern grep, [*.nc, *VR.h5]")
   
   parser.add_option("-w", "--write", dest="write",   default=False, \
                           help = "Boolean flag to write DART ascii file", action="store_true")
                           
   parser.add_option("-o", "--out",      dest="out_dir",  default="ref_files",  type="string", \
                           help = "Directory to place output files in")
                           
   parser.add_option("-p", "--plot",      dest="plot",      default=-1,  type="int",      \
                     help = "Specify a number between 0 and 20 to plot reflectivity")
 
    parser.add_option(     "--loc",      dest="loc",      default=None,  type="int",      \
                     help = "Specify location of NEWSe grid center (lat,lon)")
                  
                 
   parser.add_option(      "--thin",      dest="thin",      default=1,  type="int",      \
                     help = "Specify a number between 2 and 5 thin reflectivity")
                  
   (options, args) = parser.parse_args()

#-------------------------------------------------------------------------------

   if options.dir == None:
          
      print "\n\n ***** USER MUST SPECIFY A DIRECTORY WHERE FILES ARE *****"
      print "\n                         EXITING!\n\n"
      parser.print_help()
      print
      sys.exit(1)

   if options.thin > 1:
       _grid_dict['thin_grid'] = options.thin
      
   if options.plot < 0:
       plot_grid = False
   else:
       sweep_num = options.plot
       plot_grid = True
       plot_filename = None

   if options.realtime != None:
       year   = int(options.realtime[0:4])
       mon    = int(options.realtime[4:6])
       day    = int(options.realtime[6:8])
       hour   = int(options.realtime[8:10])
       minute = int(options.realtime[10:12])
       a_time = DT.datetime(year, mon, day, hour, minute, 0)

   g_coords = (0.5*(36.5+45.5), 0.5*(-99.5 + -91.0))

#-------------------------------------------------------------------------------

# if realtime, now find the file created within T+5 of analysis

   if options.realtime != None:

       in_filenames = get_dir_files(options.dir, options.grep, Quiet=False)

       for n, file in enumerate(in_filenames):
           f_time = DT.datetime.strptime(os.path.basename(file)[-25:-10], "%Y%m%d-%H%M%S")
           if (f_time > a_time - DT.timedelta(minutes = 2)) and (f_time <= a_time + DT.timedelta(minutes = 1)):
                last_file    = in_filenames[n]

       try:
           in_filenames = [last_file]
           print("\n prep_mrms:  RealTime FLAG is true, only processing %s\n" % (in_filenames[:]))
           rlt_filename = "%s_%s" % ("obs_seq_RF", a_time.strftime("%Y%m%d%H%M"))
       except:
           print("\n============================================================================")
           print("\n Prep_MRMS cannot find a RF file between [-2,+1] min of %s, exiting" % a_time.strftime("%Y%m%d%H%M"))
           print("\n============================================================================")
           sys.exit(1)

   else:
       in_filenames = get_dir_files(options.dir, options.grep, Quiet=False)
       out_filenames = []
       print("\n prep_mrms:  Processing %d files in the directory:  %s\n" % (len(in_filenames), options.dir))
       print("\n prep_mrms:  First file is %s\n" % (in_filenames[0]))
       print("\n prep_mrms:  Last  file is %s\n" % (in_filenames[-1]))
   
 # Make sure there is a directory to write files into....
 
   if not os.path.exists(options.out_dir):
       try:
           os.mkdir(options.out_dir)
       except:
           print("\n**********************   FATAL ERROR!!  ************************************")
           print("\n PREP_GRID3D:  Cannot create output dir:  %s\n" % options.out_dir)
           print("\n**********************   FATAL ERROR!!  ************************************")      
#-------------------------------------------------------------------------------

   for n, file in enumerate(in_filenames):
   
      print(" ================================================================================")
      print(" Constructing.....%s" % file)
             
      ref_obj = assemble_3D_grid(options.dir, os.path.basename(file), loc=g_coords, Quiet=False)

      ref_obj = dbz_masking(ref_obj, thin_zeros=_grid_dict['thin_zeros'])
      
      fsuffix = "MRMS_%s_%2.2d_KM" % (ref_obj.time.strftime('%Y%m%d%H%M'), int(ref_obj.zg[sweep_num]/100.))
      plot_filename = os.path.join(options.out_dir, fsuffix)
      grid_plot(ref_obj, sweep_num, plot_filename = plot_filename)

      fsuffix = "MRMS_%s_Composite" % (ref_obj.time.strftime('%Y%m%d%H%M'))
      plot_filename = os.path.join(options.out_dir, fsuffix)
      grid_plot(ref_obj, -1, plot_filename = plot_filename)
    
#-------------------------------------------------------------------------------
# Main program for testing...
#
if __name__ == "__main__":
    sys.exit(main())
