import time
import logging
import os
import datetime

_MRMS_directory = "/work/john.krause/realtime/grid/output"
_MRMS_obs_seq   = "/work/john.krause/realtime"

year = 2017
mon  = 03
day  = 28
hour = 18
dt        = datetime.datetime(year, mon, day, hour, 0, 0)
stop_time = datetime.datetime(year, mon, day, 21, 30, 0)
dtime     = datetime.timedelta(minutes=15)

while dt < stop_time:

     cmd = "prep_grid3d.py -d %s -w -o %s --realtime %s" % (_MRMS_directory, _MRMS_obs_seq, dt.strftime("%Y%m%d%H%M"))
     print('Running the job %s command: %s' % (time.strftime("%Y-%m-%d %H:%M:%S"), cmd))
     os.system(cmd)
     dt = dt + dtime

