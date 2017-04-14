import time
import logging
import os, sys
import datetime

_MRMS_input_directory = "/work/john.krause/realtime/grid/output"
_MRMS_obs_seq         = "./"

year       = 2017
mon        = 04
day        = [14, 14]
hour       = [20, 20]
min        = [20,  45]

plot_level = 4

start_time = datetime.datetime(year, mon, day[0], hour[0], min[0], 0)
stop_time  = datetime.datetime(year, mon, day[1], hour[1], min[1], 05)
dtime      = datetime.timedelta(minutes=15)

obs_seq_out_dir = os.path.join(_MRMS_obs_seq, start_time.strftime("%Y%m%d"))

while start_time < stop_time:

    print("\n >>>>=======BEGIN===============================================================")
    cmd = "prep_grid3d.py -d %s -w -o %s --realtime %s -p %d" % (_MRMS_input_directory, obs_seq_out_dir, start_time.strftime("%Y%m%d%H%M"), plot_level)
    print("\n prep_grid called at %s" % (time.strftime("%Y-%m-%d %H:%M:%S")))
    print(" Cmd: %s" % (cmd))
    ret = os.system("%s >> log_prep_grid3d" % cmd)
    if ret != 0:
        print("\n ============================================================================")
        print("\n Prep_Grid3D cannot find a RF file between [-2,+5] min of %s" % start_time.strftime("%Y%m%d%H%M"))
        print("\n ============================================================================")
    print("\n <<<<<=======END================================================================")

    start_time = start_time + dtime

