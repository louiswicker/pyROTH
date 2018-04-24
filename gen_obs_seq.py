#!/usr/bin/env python

import datetime as DT
import sys, os

start   = DT.datetime(2016,5,24,18,0,0)
delta_t = DT.timedelta(seconds=150)
finish  = DT.datetime(2016,5,25,3,1,0)

RunIt = True

infile = "obs_seq_2016_05_24_VR.h5"
outdir = "obs_seq_vr_files"
prefix = "obs_seq_VR"
pyDart_exe = "./pyDart.py"

if not os.path.exists(outdir):
   os.mkdir(outdir)

# copy the file to make sure we dont overwrite it

cmd = "cp %s tmp.h5" % infile
ret = os.system(cmd)

if ret != 0:
    print(" \n ERROR ERROR - temporary file could not be created - check name of input file, exiting\n")
    sys.exit(-1)

while start < finish:
    begin = start - delta_t
    end   = start + delta_t
    cmd = "%s -f tmp.h5 --start %s  --end %s --search --hdf2ascii" % (pyDart_exe, 
                                                                  begin.strftime("%Y,%m,%d,%H,%M,%S"), \
                                                                  end.strftime("%Y,%m,%d,%H,%M,%S"))
    print(cmd)
    if RunIt:
       os.system(cmd)

    cmd = "mv tmp.tmp.out %s/%s_%s.out" % (outdir,prefix, start.strftime("%Y%m%d_%H%M%S"))

    print(cmd)
    if RunIt:
       os.system(cmd)

    start = start + DT.timedelta(seconds=900)


