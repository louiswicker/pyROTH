#!/usr/bin/env python

from apscheduler.schedulers.blocking import BlockingScheduler
import time
import logging
import os,sys
import datetime
from MRMS_2_Dart import run_Prep_Grid3d, run_MRMS, kill_MRMS_Programs

today = time.strftime("%Y%m%d")

start = datetime.datetime(2017,3,31,11,55,0)
end   = datetime.datetime(2017,3,31,13,06,0)

print("\n==============================================================================")
print("\n Starting cron_MRMS_Dart script at: %s " % time.strftime("%Y-%m-%d %H:%M:%S"))
print("\n Todays output will be written into the directory: %s " % today)
print("==============================================================================\n")

logging.basicConfig()
sched = BlockingScheduler()

# schedule the main mrms programs
@sched.scheduled_job('cron', hour=11, minute=50)
def scheduled_job():
   print("\n  ----> Begin run MRMS processing at %s <----- \n" % time.strftime("%Y:%m:%d %H:%M"))

@sched.scheduled_job('cron', hour=13, minute=15)
def scheduled_job():
   print("\n  ----> End processing at %s <----- \n" % time.strftime("%Y:%m:%d %H:%M"))

# schedule the prep_grid3d...

@sched.scheduled_job('cron', start_date=start, end_date=end, minute="5,20,35,50")
def scheduled_job():
   print("\n  ----> Running prep_grid3d.py at %s <----- \n" % time.strftime("%Y:%m:%d %H:%M"))

sched.start()
