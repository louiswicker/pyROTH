# Author: Jonathan J. Helmus (jhelmus@anl.gov)
# License: BSD 3 clause

import scipy
import copy
import matplotlib.pyplot as plt
import pyart
import matplotlib
import pylab as P
from scipy import signal
import numpy as N
import math
import sys
import netCDF4
from optparse import OptionParser
from netCDF4 import num2date
from netcdftime import utime
import os
import ctables
from mpl_toolkits.basemap import Basemap
from pyproj import Proj
import time as timeit
import numpy.ma as ma
import netCDF4 as ncdf
import glob
#from Plotting import shapefiles
shapefiles = None

#---------------------------------------------------------------------------------------------------


# Lambert conformal stuff
tlat1 = 35.3500
tlat2 = 35.8900
cen_lat = 35.6200
cen_lon = -97.9120

_counties = True
_states   = False

# Colorscale information
ref_scale = (0.,74.)
vr_scale  = (-40.,40.)

# Range rings in km
range_rings = [25, 50, 75, 100, 125] 

x_plot_range = [-150000., 150000.]
y_plot_range = [-150000., 150000.] 

# Debug print stuff....
debug = False

#======================= EXAMPLE for Container CLASS ====================================
# class data_bin(object):
#     pass
#    new_obj = data_bin()
#    new_obj.fields = {}
#    new_obj.loc = copy.deepcopy(display_obj.loc)
#    new_obj.x = copy.deepcopy(display_obj.x)
#    new_obj.y = copy.deepcopy(display_obj.y)


#===============================================================================
def mybasemap(glon, glat, r_lon, r_lat, scale = 1.0, ticks = True, 
              resolution='c',area_thresh = 10., shape_env = False, 
              counties=False, states = False, lat_lines=True, lon_lines=True, 
              pickle = False, ax=None):

   tt = timeit.clock()
   
   map = Basemap(llcrnrlon=glon.min(), llcrnrlat=glat.min(), \
                  urcrnrlon=glon.max(), urcrnrlat=glat.max(), \
                  lat_0 = r_lat, lon_0=r_lon, \
                  projection = 'lcc',      \
                  resolution=resolution,   \
                  area_thresh=area_thresh, \
                  suppress_ticks=ticks, ax=ax)

   if counties or _counties:
      map.drawcounties()
      
   if states or _states:
       map.drawstates()
       
   if lat_lines == True:
       lat_lines = N.arange(30., 60., 0.5)
       map.drawparallels(lat_lines, labels=[True, False, False, False])
       
   if lon_lines == True:
       lon_lines = N.arange(-110., 70., 0.5)
       map.drawmeridians(lon_lines, labels=[False, False, False, True])
       
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
def create_ppi_map(radar, xr, yr, plot_range_rings=True, ax=None, **kwargs):

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
      for ring in range_rings:
         xpts = radar_x + ring * 1000. * N.sin(angle)
         ypts = radar_y + ring * 1000. * N.cos(angle)
         map.plot(xpts, ypts, color = 'gray', alpha = 0.5, linewidth = 1.0, ax=ax)  
                 
   return map, xmap, ymap, radar_x, radar_y
   
#############################################################################################
def plot_ppi_map(radar, field, level= 0, cmap=ctables.Carbone42, vRange=None, var_label = None, \
                 plot_range_rings=True, ax=None, **kwargs):
                 
   start = radar.get_start(level)
   end   = radar.get_end(level) + 1
   data  = radar.fields[field]['data'][start:end]
   
# Fix super-res issue?

   if N.sum(data.mask == True) == data.size:
      start = radar.get_start(level+1)
      end   = radar.get_end(level+1) + 1
      data  = radar.fields[field]['data'][start:end]

   xr    = radar.gate_x['data'][start:end] 
   yr    = radar.gate_y['data'][start:end] 
   
   map, x, y, radar_x, radar_y = create_ppi_map(radar, xr, yr, plot_range_rings=True, ax=ax, **kwargs)
   
   plot = map.pcolormesh(x, y, data, cmap=cmap, vmin = vRange[0], vmax = vRange[1], ax=ax)
   cbar = map.colorbar(plot, location = 'right', pad = '3%')
   
   ax.set_xlim(radar_x+x_plot_range[0],radar_x+x_plot_range[1])
   ax.set_ylim(radar_y+y_plot_range[0],radar_y+y_plot_range[1])
   
   if var_label:
       cbar.set_label(var_label, fontweight='bold')
   else:
       cbar.set_label(field, fontweight='bold')

   time = radar.time['units'].split(" ")[2]
   time = time.replace("-","_")
   time = time.replace("T"," ")
   el   = radar.elevation['data'][level]
   
   title_string = "%s    %s    EL:  %4.2f deg" % (field, time, el)

   P.title(title_string)
   
   print("\nCompleted plot for %s" % field)
     
   return 
     
#############################################################################################
# MAIN program

filenames = ["/Users/Louis.Wicker/Software/PythonProjects/pyROTH/KTLX/KTLX20130520_205045_V06"]

print filenames[0]
print filenames[-1]

output_dir = "testimages"

out_filenames = []

if filenames[0][-3:] == "V06":
  for item in filenames:
    strng = os.path.basename(item)[0:17]
    strng = strng[0:4] + "_" + strng[4:]
    strng = os.path.join(output_dir, strng)
    out_filenames.append(strng)

for n, filename in enumerate(filenames):

  fig, axes = P.subplots(1, 2, sharey=True, figsize=(15,6))
  
  radar = pyart.io.read_nexrad_archive(filename)

  outfile  = plot_ppi_map(radar, "reflectivity", vRange=ref_scale, cmap=ctables.NWSRef, \
                          ax=axes[0], var_label='Reflectivity', shape_env=shapefiles)

  outfile  = plot_ppi_map(radar, "velocity", vRange=(-40,40), cmap=ctables.Carbone42, ax=axes[1], 
                          var_label='Radial Velocity', shape_env=shapefiles)

  fig.subplots_adjust(left=0.06, right=0.90, top=0.90, bottom=0.1, wspace=0.35)
  
#   plot range rings at 10, 20, 30 and 40km
#   display.plot_range_rings([10, 20, 30, 40])
#   
#   plots cross hairs
#   display.plot_line_xy(np.array([-40000.0, 40000.0]), np.array([0.0, 0.0]),
#                        line_style='k-')
#   display.plot_line_xy(np.array([0.0, 0.0]), np.array([-20000.0, 200000.0]),
#                        line_style='k-')
# 
# Indicate the radar location with a point
#   display.plot_point(radar.longitude['data'][0], radar.latitude['data'][0])

  print("\n Saving file:  %s.png" % (out_filenames[n]))

  P.savefig("%s.png" % out_filenames[n], format="png", dpi=300)

  P.show()
# 
# 
# # 		outfile  = plot_ppi_map(display_obj, "reflectivity", vRange=ref_scale, cmap=ctables.NWSRef, \
# # 		                        ax=axes[0], var_label='Reflectivity', shape_env=shapefiles)
# # 
# # 		outfile  = plot_ppi_map(display_obj, "velocity", vRange=(-40,40), cmap=ctables.Carbone42, ax=axes[1], 
# # 		                        var_label='Radial Velocity', shape_env=shapefiles)
# # 
# # 		
# 
#P.show()

