#!/usr/bin/env python

from apscheduler.schedulers.blocking import BlockingScheduler
import time
import os,sys
import datetime
from MRMS_2_Dart import run_Prep_Grid3d, run_MRMS, kill_MRMS_Programs

today = time.strftime("%Y%m%d")

print("\n==============================================================================")
print("\n Starting cron_MRMS_Dart script at: %s " % time.strftime("%Y-%m-%d %H:%M:%S"))
print("\n Todays output will be written into the directory: %s " % today)
print("==============================================================================\n")

sched = BlockingScheduler()

# schedule the prep_grid3d...

@sched.scheduled_job('cron', minute="6,20,35,50")
def scheduled_job():
   print("\n  ----> Running prep_grid3d.py at %s <----- \n" % time.strftime("%Y:%m:%d %H:%M"))
   today = time.strftime("%Y%m%d")
   ret = run_Prep_Grid3d(today)
   if ret != None:
      print("\n  ----> prep_grid3d did not return None:  %s <----- \n" % ret)

sched.start()
