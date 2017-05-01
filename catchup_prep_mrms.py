#!/usr/bin/env python


import time
import logging
import os, sys
import datetime

_MRMS_feed       = "/work/LDM/MRMS/"
_MRMS_obs_seq    = "/work/wicker/REALTIME/"
_NEWSe_grid_info = "/scratch/wof/realtime/radar_files"

year       = 2017
mon        = 05
day        = [1, 1]
hour       = [18, 19]
min        = [30,  30]

plot_level = 3

start_time = datetime.datetime(year, mon, day[0], hour[0], min[0], 0)
stop_time  = datetime.datetime(year, mon, day[1], hour[1], min[1], 05)
dtime      = datetime.timedelta(minutes=15)

obs_seq_out_dir = os.path.join(_MRMS_obs_seq, start_time.strftime("%Y%m%d"))

# create path for NEWSe radar file

radar_csh_file = os.path.join(_NEWSe_grid_info, ("radars.%s.csh" % start_time.strftime("%Y%m%d")))
print radar_csh_file

# Parse center lat and lon out of the c-shell radar file - HARDCODED!
# If the file does not exist, then we exit out of this run

try:
    fhandle = open(radar_csh_file)
except:
    print("\n ============================================================================")
    print("\n CANNOT OPEN radar CSH file, exiting MRMS processing:  %s" % radar_csh_file)
    print("\n ============================================================================")
    sys.exit(1)

all_lines  = fhandle.readlines()
lat = float(all_lines[7].split(" ")[2])
lon = float(all_lines[8].split(" ")[2])
fhandle.close()

print("\n ============================================================================")
print("\n Lat: %f  Lon: %f centers will be used for MRMS sub-domain" % (lat, lon))
print("\n ============================================================================")

# Main catchup loop

while start_time < stop_time:

    MRMS_dir = os.path.join(_MRMS_feed, start_time.strftime("%Y/%m/%d"))
    
    print("\n Reading from operational MRMS directory:  %s\n" % MRMS_dir)
    
    print("\n >>>>=======BEGIN===============================================================")
    cmd = "prep_mrms.py -d %s -w -o %s --realtime %s -p %d --loc %f %f"  %  \
          (MRMS_dir, obs_seq_out_dir, start_time.strftime("%Y%m%d%H%M"), plot_level, lat, lon)

    print("\n Prep_MRMS called at %s" % (time.strftime("%Y-%m-%d %H:%M:%S")))
    print(" Cmd: %s" % (cmd))
    ret = os.system("%s >> log_prep_MRMS" % cmd)
    if ret != 0:
        print("\n ============================================================================")
        print("\n Prep_MRMS cannot find a RF file between [-2,+1] min of %s" % start_time.strftime("%Y%m%d%H%M"))
        print("\n ============================================================================")
    print("\n <<<<<=======END================================================================")

    start_time = start_time + dtime

