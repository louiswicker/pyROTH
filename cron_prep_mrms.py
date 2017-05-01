#!/usr/bin/env python

from apscheduler.schedulers.blocking import BlockingScheduler
import time
import os
import sys
import datetime

_MRMS_feed       = "/work/LDM/MRMS/"
_MRMS_obs_seq    = "/work/wicker/REALTIME/"
_NEWSe_grid_info = "/scratch/wof/realtime/radar_files"
_NEWSe_prep_mrms = "/work/wicker/REALTIME/pyroth/prep_mrms.py"

plot_level = 3

#-----------------------------------------------------------------------------
# Utility to round the datetime object to nearest 15 min....
def quarter_datetime(dt):
    minute = (dt.minute//15)*15
    return dt.replace(minute=0, second=0)+datetime.timedelta(minutes=minute)

# Main program
# Okay, some fudge here:  need to keep the today as today according to local time..

local_today = time.localtime()
today = "%s%2.2d%2.2d" % (local_today.tm_year, local_today.tm_mon, local_today.tm_mday)

print("\n ==============================================================================")
print("\n Starting cron_prep_MRMS script at: %s " % time.strftime("%Y-%m-%d %H:%M:%S"))
print("\n Todays output will be written into the directory: %s " % today)
print(" ==============================================================================\n")

obs_seq_out_dir = os.path.join(_MRMS_obs_seq, today)

# create path for NEWSe radar file

radar_csh_file = os.path.join(_NEWSe_grid_info, ("radars.%s.csh" % today))

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

sched = BlockingScheduler()

# schedule the prep_grid3d...

@sched.scheduled_job('cron', minute="5,20,35,50")
def scheduled_job():

    gmt = time.gmtime()  # for file names, here we need to use GMT time

    dt = quarter_datetime(datetime.datetime(*gmt[:6]))

    str_time = dt.strftime("%Y%m%d%H%M")

    MRMS_dir = os.path.join(_MRMS_feed, dt.strftime("%Y/%m/%d"))
    
    print("\n Reading from operational MRMS directory:  %s\n" % MRMS_dir)
    
    print("\n >>>>=======BEGIN===============================================================")
    cmd = "%s -d %s -w -o %s --realtime %s -p %d --loc %f %f"  %  \
          (_NEWSe_prep_mrms, MRMS_dir, obs_seq_out_dir, dt.strftime("%Y%m%d%H%M"), plot_level, lat, lon)

    print("\n Prep_MRMS called at %s" % (dt.strftime("%Y-%m-%d %H:%M:%S")))
    print(" Cmd: %s" % (cmd))

    ret = os.system("%s >> log_prep_MRMS" % cmd)
    if ret != 0:
        print("\n ============================================================================")
        print("\n Prep_MRMS cannot find a RF file between [-2,+1] min of %s" % dt_time.strftime("%Y%m%d%H%M"))
        print("\n ============================================================================")

    print("\n <<<<<=======END================================================================")

sched.start()
