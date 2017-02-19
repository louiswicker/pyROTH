#!/usr/bin/env python

import datetime as DT
import sys, os

start   = DT.datetime(2016,5,9,18,0,0)
delta_t = DT.timedelta(seconds=150)
finish  = DT.datetime(2016,5,10,3,1,0)

RunIt = True

infile = "../09May/VR_20DBZ/obs_seq_2016_05_09_VR.h5"
oufile = "vr_files/obs_seq_VR"
pyDart_exe = "./pyDart.py"

# copy the file to make sure we dont overwrite it

cmd = "cp %s tmp.h5" % infile

while start < finish:
    begin = start - delta_t
    end   = start + delta_t
    cmd = "%s -f tmp.h5 --start %s  --end %s --search --hdf2ascii" % (pyDart_exe, 
                                                                  begin.strftime("%Y,%m,%d,%H,%M,%S"), \
                                                                  end.strftime("%Y,%m,%d,%H,%M,%S"))
    print(cmd)
    if RunIt:
       os.system(cmd)

    cmd = "mv tmp.tmp.out %s_%s.out" % (oufile, start.strftime("%Y%m%d_%H%M%S"))

    print(cmd)
    if RunIt:
       os.system(cmd)

    start = start + DT.timedelta(seconds=900)


