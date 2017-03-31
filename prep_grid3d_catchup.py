import time
import logging
import os
import datetime

_MRMS_input_directory = "/work/john.krause/realtime/grid/output"
_MRMS_obs_seq   = "/work/john.krause/realtime"

year = 2017
mon  = 03
day  = 31
hour = 03
dt        = datetime.datetime(year, mon, day, hour, 0, 0)
stop_time = datetime.datetime(year, mon, day, 03, 31, 0)
dtime     = datetime.timedelta(minutes=15)

day       = datetime.datetime(year, mon, day, 0, 0)

obs_seq_out_dir = os.path.join(_MRMS_obs_seq, time.strftime("%Y%m%d"))

while dt < stop_time:

    cmd = "prep_grid3d.py -d %s -w -o %s --realtime %s -p 4" % (_MRMS_input_directory, obs_seq_out_dir, dt.strftime("%Y%m%d%H%M"))
    print("\n Run_MRMS running job: at %s" % (time.strftime("%Y-%m-%d %H:%M:%S")))
    print(" Command: %s\n" % (cmd))
    os.system("%s >> log_prep_grid3d" % cmd)

    dt = dt + dtime

