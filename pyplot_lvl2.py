#!/usr/bin/env python

import matplotlib.pyplot as plt
import pyart
import matplotlib
import pylab as P
import numpy as N
import numpy as np
import sys
import glob
from optparse import OptionParser
import os
import ctables
import time as timeit

from mpl_toolkits.basemap import Basemap
from pyproj import Proj

#---------------------------------------------------------------------------------------------------
#
# Useful defaults if you are doing a lot of plot and dont want to put them on the command line
#
_counties = True
_states   = False

_shapefiles = None
_level      = 0

_min_dbz    = 35.

# Colorscale information
_ref_scale = (0.,74.)
_vr_scale  = (-40.,40.)

# Colortables
_ref_ctable = ctables.NWSRef
_vr_ctable  = ctables.Carbone42

# Range rings in km
_plot_RangeRings = False
_range_rings = [25, 50, 75, 100, 125] 

# How much to plot
_max_range    = 150000.
_x_plot_range = [-_max_range, _max_range]
_y_plot_range = [-_max_range, _max_range] 
_lat_lines    = True
_lon_lines    = True

# Output file format
_file_type = "png"

# this turns on some timing info
debug = False

########################################################################
# Volume Prep:  Threshold data and unfold velocities

def volume_prep(radar, unfold_type="phase"):
          
# Mask data beyond max_range

  max_range_gate = np.abs(radar.range['data'] - _max_range).argmin()
  radar.fields['reflectivity']['data'][:,max_range_gate:] = np.ma.masked
  radar.fields['velocity']['data'][:,max_range_gate:] = np.ma.masked

# Filter based on masking, dBZ threshold, and invalid gates

  gatefilter = pyart.correct.GateFilter(radar)
  gatefilter.exclude_invalid('velocity')
  gatefilter.exclude_invalid('reflectivity')
  gatefilter.exclude_masked('velocity')
  gatefilter.exclude_masked('reflectivity')
  gatefilter.exclude_below('reflectivity', _min_dbz)

  pyart.correct.despeckle.despeckle_field(radar, 'velocity', threshold=-100, size=10, gatefilter=gatefilter, delta=5.0)
  pyart.correct.despeckle.despeckle_field(radar, 'reflectivity', threshold=-100, size=10, gatefilter=gatefilter, delta=5.0)

# Dealias the velocity data

  if unfold_type == "phase":
    dealiased_radar = pyart.correct.dealias_unwrap_phase(radar, unwrap_unit='sweep', 
                                   nyquist_vel=None, check_nyquist_uniform=True, 
                                   gatefilter=gatefilter, rays_wrap_around=None, 
                                   keep_original=False, set_limits=True, 
                                   vel_field='velocity', corr_vel_field=None, 
                                   skip_checks=False)
                                   
    radar.add_field('unfolded velocity', dealiased_radar)
    
    return True
                                   

  if unfold_type == "region":
    dealiased_radar = pyart.correct.dealias_region_based(radar, interval_splits=3, 
                                   interval_limits=None, skip_between_rays=100, 
                                   skip_along_ray=100, centered=True, 
                                   nyquist_vel=None, check_nyquist_uniform=True, 
                                   gatefilter=gatefilter, rays_wrap_around=None, 
                                   keep_original=False, set_limits=True, 
                                   vel_field='velocity', corr_vel_field=None)
     
    radar.add_field('unfolded velocity', dealiased_radar)

    return True
    
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

# If it gets here, there was a problem

  return False
#===============================================================================
def mybasemap(glon, glat, r_lon, r_lat, scale = 1.0, supress_ticks = True, 
              resolution='c',area_thresh = 10., shape_env = False, 
              counties=False, states = False, lat_lines=_lat_lines, lon_lines=_lon_lines, 
              pickle = False, ax=None, **kwargs):

   tt = timeit.clock()
   
   map = Basemap(llcrnrlon=glon.min(), llcrnrlat=glat.min(), \
                  urcrnrlon=glon.max(), urcrnrlat=glat.max(), \
                  lat_0 = r_lat, lon_0=r_lon, \
                  projection = 'lcc',      \
                  resolution=resolution,   \
                  area_thresh=area_thresh, \
                  suppress_ticks=supress_ticks, ax=ax)

   if counties or _counties:
      map.drawcounties()
      
   if states or _states:
       map.drawstates()
       
   if lat_lines == True:
       lat_lines = N.arange(30., 60., 0.5)
       map.drawparallels(lat_lines, labels=[1,0,0,0], fontsize=10)
       
   if lon_lines == True:
       lon_lines = N.arange(-110., 70., 0.5)
       map.drawmeridians(lon_lines, labels=[0,0,0,1], fontsize=10)
       
   # Shape file stuff

   if shape_env:

      for item in shape_env:
          items = item.split(",")
          shapefile  = items[0]
          color      = items[1]
          linewidth  = float(items[2])

          s = map.readshapefile(shapefile,'shapeinfo',drawbounds=False)

          for shape in map.shapeinfo:
              xx, yy = zip(*shape)
              map.plot(xx,yy,color=color,linewidth=linewidth,ax=ax)

   # pickle the class instance.

   if debug:  print(timeit.clock()-tt,' secs to create original Basemap instance')

   if pickle:
      pickle.dump(map,open('mymap.pickle','wb'),-1)
      print(timeit.clock()-tt,' secs to create original Basemap instance and pickle it')

   return map

#############################################################################################
def create_ppi_map(radar, xr, yr, plot_range_rings=_plot_RangeRings, ax=None, **kwargs):

# Lambert conformal defaults
   tlat1 = 35.3500
   tlat2 = 35.8900

   radar_lon = radar.longitude['data'][0]
   radar_lat = radar.latitude['data'][0]

   p1 = Proj(proj='lcc', ellps='WGS84', datum='WGS84', lat_1=tlat1, lat_2=tlat2, lon_0=radar_lon)

   x_offset, y_offset = p1(radar_lon, radar_lat)
   
   x = xr + x_offset
   y = yr + y_offset

   lon, lat = p1(x, y, inverse=True)

   map = mybasemap(lon, lat, radar_lon, radar_lat, ax=ax, **kwargs)

   radar_x, radar_y = map(radar_lon, radar_lat)
   xmap, ymap = map(lon, lat)

   if plot_range_rings:
      angle = N.linspace(0., 2.0 * N.pi, 360)
      for ring in _range_rings:
         xpts = radar_x + ring * 1000. * N.sin(angle)
         ypts = radar_y + ring * 1000. * N.cos(angle)
         map.plot(xpts, ypts, color = 'gray', alpha = 0.5, linewidth = 1.0, ax=ax)  
                 
   return map, xmap, ymap, radar_x, radar_y
   
#############################################################################################
def plot_ppi_map(radar, field, level = 0, cmap=ctables.Carbone42, vRange=None, var_label = None, \
                 plot_range_rings=True, ax=None, zoom = None, **kwargs):
                 
   start = radar.get_start(level)
   end   = radar.get_end(level) + 1
   data  = radar.fields[field]['data'][start:end]
   
# Fix super-res blank velocity field to next scan level...

   if N.sum(data.mask == True) == data.size:
      start = radar.get_start(level+1)
      end   = radar.get_end(level+1) + 1
      data  = radar.fields[field]['data'][start:end]

   xr    = radar.gate_x['data'][start:end] 
   yr    = radar.gate_y['data'][start:end] 
   
   map, x, y, radar_x, radar_y = create_ppi_map(radar, xr, yr, ax=ax, **kwargs)
   
   plot = map.pcolormesh(x, y, data, cmap=cmap, vmin = vRange[0], vmax = vRange[1], ax=ax)
   cbar = map.colorbar(plot, location = 'right', pad = '3%')
   
   if zoom:
      ax.set_xlim(radar_x+1000.*zoom[0],radar_x+1000.*zoom[1])
      ax.set_ylim(radar_y+1000.*zoom[2],radar_y+1000.*zoom[3])
   else:   
      ax.set_xlim(radar_x+_x_plot_range[0],radar_x+_x_plot_range[1])
      ax.set_ylim(radar_y+_y_plot_range[0],radar_y+_y_plot_range[1])

   time = radar.time['units'].split(" ")[2]
#time = time.replace("-","_")
   time = time.replace("T"," ")
   el   = radar.elevation['data'][level]
   
   if var_label != None:
       cbar.set_label(var_label, fontweight='bold')
       title_string = "%s    EL:  %4.2f deg" % (var_label, el)
   else:
       cbar.set_label(field, fontweight='bold')
       title_string = "%s    EL:  %4.2f deg" % (field,  el)

   P.suptitle(time, fontsize=16)
   P.title(title_string)
   
   print("\nCompleted plot for %s" % field)
     
   return 
     
########################################################################
# Main function

if __name__ == "__main__":

  print ' ================================================================================'
  print ''
  print ''
  print '                   BEGIN PROGRAM pyPlot_lvl2                    '
  print ''

  parser = OptionParser()
  parser.add_option("-d", "--dir",       dest="dname",     default=None,  type="string", \
                    help = "Directory of files to process")
                 
  parser.add_option("-o", "--out",       dest="out_dir",     default="images",  type="string", \
                    help = "Directory to place png files in")
                 
  parser.add_option("-f", "--file",      dest="fname",     default=None,  type="string", \
                    help = "filename of NEXRAD level II volume to process")

  parser.add_option("-l", "--level",      dest="level",     default=0,  type="int", \
                    help = "Tilt index of NEXRAD sweep--> 0: 0.5 degrees")
                 
                 
  parser.add_option("-u", "--unfold",    dest="unfold",    default="phase",  type="string", \
                    help = "dealiasing method to use (phase or region, default = phase)")
  
  parser.add_option("-i", "--interactive", dest="interactive", default=False,  action="store_true",     \
                     help = "Boolean flag to plot image to screen")  
                     
  parser.add_option("-s", "--shapefiles", dest="shapefiles", default=None, type="string",    \
                     help = "Name of system env shapefile you want to add to the plots.")
                     
  parser.add_option(       "--zoom",     dest="zoom", type="int", nargs = 4, \
                     help="zoom in on (in KM) of plot 4 args required: xmin xmax ymin ymax")
                     
  (options, args) = parser.parse_args()
  
  print ''
  print ' ================================================================================'
  
  if not os.path.exists(options.out_dir):
    os.mkdir(options.out_dir)

  out_filenames = []
  in_filenames  = []

  if options.dname == None:
          
    if options.fname == None:
      print "\n\n ***** USER MUST SPECIFY NEXRAD LEVEL II (MESSAGE 31) FILE OR DIRECTORY! *****"
      print "\n                         EXITING!\n\n"
      parser.print_help()
      print
      sys.exit(1)
      
    else:
      in_filenames.append(os.path.abspath(options.fname))
      strng = os.path.basename(in_filenames[0])[0:-3]
      strng = strng[0:4] + "_" + strng[4:]
      strng = os.path.join(options.out_dir, strng)
      out_filenames.append(strng) 

  else:
    in_filenames = glob.glob("%s/*" % os.path.abspath(options.dname))
    print("\n pyROTH:  Processing %d files in the directory:  %s\n" % (len(in_filenames), options.dname))
    print("\n pyROTH:  First file is %s\n" % (in_filenames[0]))
    print("\n pyROTH:  Last  file is %s\n" % (in_filenames[-1]))

    if in_filenames[0][-3:] == "V06":
      for item in in_filenames:
        strng = os.path.basename(item)[0:17]
        strng = strng[0:4] + "_" + strng[4:]
        strng = os.path.join(output_dir, strng)
        out_filenames.append(strng)

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
      
  if options.shapefiles:
      shapefiles = options.shapefiles
  else:
      shapefiles = _shapefiles

#===========================
# Main processing loop....
      
  for n, filename in enumerate(in_filenames):

    fig, axes = P.subplots(1, 2, sharey=True, figsize=(15,6))
  
    volume = pyart.io.read_nexrad_archive(filename)
    
  # unfolding can fail if the data are not written quite write - instead of quiting, try to unfold with region method
    if unfold_type == None:
        vr_field = "velocity"
        vr_label = "Radial Velocity"

    else:
        try:
            gatefilter = volume_prep(volume, unfold_type=unfold_type) 
            vr_field = "unfolded velocity"
            vr_label = "Unfolded Radial Velocity"
        except:
            try:
                print("\n ----> Phase unfolding method has failed!! Trying region unfolding method\n")
                gatefilter = volume_prep(volume, unfold_type="region")
                vr_field = "unfolded velocity"
                vr_label = "Unfolded Radial Velocity"
            except:
                print("\n ----> Both unfolding methods have failed!! Turning off unfolding\n\n")
                vr_field = "velocity"
                vr_label = "Radial Velocity"

    outfile  = plot_ppi_map(volume, "reflectivity", level=options.level, vRange=_ref_scale, cmap=_ref_ctable, \
                            ax=axes[0], var_label='Reflectivity', shape_env=shapefiles, zoom=options.zoom)

    outfile  = plot_ppi_map(volume, vr_field, level=options.level, vRange=_vr_scale, cmap=_vr_ctable, ax=axes[1], 
                            var_label=vr_label, shape_env=shapefiles, zoom=options.zoom)

    fig.subplots_adjust(left=0.06, right=0.90, top=0.90, bottom=0.1, wspace=0.35)
    
    outfile = "%s%4.2fDEG" % (out_filenames[n], volume.elevation['data'][options.level])
  
    print("\n Saving file:  %s.png" % (outfile))

    P.savefig("%s.%s" % (outfile, _file_type), format=_file_type, dpi=300)

    if options.interactive:
        P.show()
