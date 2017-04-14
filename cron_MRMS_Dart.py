#!/usr/bin/env python

from apscheduler.schedulers.blocking import BlockingScheduler
import time
import logging
import os,sys
import datetime
from MRMS_2_Dart import run_Prep_Grid3d, run_MRMS, kill_MRMS_Programs

time_struct = time.localtime()
year  = time_struct.tm_year
month = time_struct.tm_mon
day   = time_struct.tm_mday
mon   = 4
day   = 1

start = datetime.datetime(year, month, day, 9, 1, 0)
end   = datetime.datetime(year, month, day ,21, 06, 0)

print("\n==============================================================================")
print("\n Starting cron_MRMS_Dart script at: %s " % time.strftime("%Y-%m-%d %H:%M:%S"))
print("==============================================================================\n")

logging.basicConfig()
sched = BlockingScheduler()

# schedule the main mrms programs

@sched.scheduled_job('cron', start_date=start, end_date=end, hour=11, minute=12)
def scheduled_job():
   print("\n  ----> Begin run MRMS processing at %s <----- \n" % time.strftime("%Y:%m:%d %H:%M"))
   ret = run_MRMS()

@sched.scheduled_job('cron', start_date=start, end_date=end, hour=21, minute=05)
def scheduled_job():
   print("\n  ----> End processing at %s <----- \n" % time.strftime("%Y:%m:%d %H:%M"))
   ret = kill_MRMS_Programs()


# schedule the prep_grid3d...

@sched.scheduled_job('cron', start_date=start, end_date=end, minute="5,20,35,50")
def scheduled_job():
   print("\n  ----> Running prep_grid3d.py at %s <----- \n" % time.strftime("%Y:%m:%d %H:%M"))
   today = time.strftime("%Y%m%d")
   run_Prep_Grid3d(today)

# Get started!

sched.start()
