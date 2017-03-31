import time
import logging
import os
import datetime

_MRMS_input_directory = "/work/john.krause/realtime/grid/output"
_MRMS_obs_seq   = "/work/john.krause/realtime"

year = 2017
mon  = 03
day  = 31
hour = 02
dt        = datetime.datetime(year, mon, day, hour, 0, 0)
stop_time = datetime.datetime(year, mon, day, 03, 31, 0)
dtime     = datetime.timedelta(minutes=15)

day       = datetime.datetime(year, mon, day, 0, 0)

obs_seq_out_dir = os.path.join(_MRMS_obs_seq, time.strftime("%Y%m%d"))

while dt < stop_time:

    print("\n >>>>=======BEGIN===============================================================")
    cmd = "prep_grid3d.py -d %s -w -o %s --realtime %s -p 4" % (_MRMS_input_directory, obs_seq_out_dir, dt.strftime("%Y%m%d%H%M"))
    print("\n prep_grid called at %s" % (time.strftime("%Y-%m-%d %H:%M:%S")))
    print(" Cmd: %s" % (cmd))
    ret = os.system("%s >> log_prep_grid3d" % cmd)
    if ret != 0:
        print("\n ============================================================================")
        print("\n Prep_Grid3D cannot find a RF file between [-1,+6] min of %s" % dt.strftime("%Y%m%d%H%M"))
        print("\n ============================================================================")
    print("\n <<<<<=======END================================================================")

    dt = dt + dtime

