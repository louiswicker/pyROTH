#!/usr/bin/env python

from apscheduler.schedulers.blocking import BlockingScheduler
import time
import logging
import os,sys
import datetime
from MRMS_2_Dart import run_Prep_Grid3d, run_MRMS, kill_MRMS_Programs

today = time.strftime("%Y%m%d")

print("\n==============================================================================")
print("\n Starting cron_MRMS_Dart script at: %s " % time.strftime("%Y-%m-%d %H:%M:%S"))
print("\n Todays output will be written into the directory: %s " % today)
print("==============================================================================\n")

logging.basicConfig()
sched = BlockingScheduler()

# schedule the main mrms programs
#you need to move output and directories around, right?

# schedule the prep_grid3d...

@sched.scheduled_job('cron', minute="5,20,35,50")
def scheduled_job():
    run_Prep_Grid3d(today)

sched.start()
