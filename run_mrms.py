from apscheduler.schedulers.blocking import BlockingScheduler
import time as timeit
import logging
import os
_MRMS_directory = "/work/john.krause/realtime/grid/output"

logging.basicConfig()

sched = BlockingScheduler()

@sched.scheduled_job('interval', seconds=60)
def timed_job():
    print('This job is run every 60 seconds.')

@sched.scheduled_job('cron', minute="5,20,35,50")
def scheduled_job():
    cmd = "prep_grid3d.py -d %s -l -w" % _MRMS_directory
    print('Running the job %s command: %s' % (timeit.strftime("%Y-%m-%d %H:%M:%S"), cmd))
#   os.system(cmd)

sched.start()
