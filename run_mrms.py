#!/usr/bin/env python

from apscheduler.schedulers.blocking import BlockingScheduler
import time
import logging
import os
import datetime

_MRMS_directory = "/work/john.krause/realtime/grid/output"
_MRMS_obs_seq   = "/work/john.krause/realtime"

def quarter_datetime(dt):
    minute = (dt.minute//15)*15
    return dt.replace(minute=0, second=0)+datetime.timedelta(minutes=minute)

logging.basicConfig()

sched = BlockingScheduler()

@sched.scheduled_job('cron', minute="5,20,35,50")
def scheduled_job():
    gmt = time.gmtime()
    dt = quarter_datetime(datetime.datetime(*gmt[:6]))
    str_time = dt.strftime("%Y%m%d%H%M")
    cmd = "prep_grid3d.py -d %s -w -o %s --realtime %s" % (_MRMS_directory, _MRMS_obs_seq, str_time)
    print("\n Run_MRMS running job: at %s" % (time.strftime("%Y-%m-%d %H:%M:%S")))
    print(" Command: %s\n" % (cmd))
    os.system("%s >> log_prep_grid3d" % cmd)

sched.start()
