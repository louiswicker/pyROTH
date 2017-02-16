#!/usr/bin/env python
#############################################################
#
#  script to download large numbers of NEXRAD LVL2 files 
#  from Amazon S3 and process them to create VR superob files 
#  in HDF5 and DART obs_seq files
#
#  Extra Python requirements:  a boto3 install
#       this is the software interface to AWS
#
#############################################################
#
# Written by Lou Wicker Feb, 2017
#
# Thanks to Tony Reinhart for getting me the code for S3
#
#############################################################

import datetime as DT
import sys, os
import boto3
from optparse import OptionParser
import datetime as DT
from multiprocessing import Pool
import time as cpu

#=======================================================================================================================
# Definitions

_nthreads = 2

_chk_dir_size  = 300

_region        = 'us-east-1'

_delta_t       = DT.timedelta(hours=1)  

_wget_string   = "wget -c https://noaa-nexrad-level2.s3.amazonaws.com/"

_prep_string   = "prep_nexrad.py --start %s --end %s %s -r %s >& out_prep_%s %s "

debug = True

#=======================================================================================================================
# RunProcess is a function that runs a system command for parallel processing

def RunProcess(cmd):

    print("\n Executing command:  %s " % cmd)

    os.system(cmd)

    print("\n %s completed...." % cmd)

    return

#=======================================================================================================================
def get_folder_size(start_path = '.'):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    return total_size

#=======================================================================================================================
# Parse and run NEWS csh radar file

def parse_NEWSe_radar_file(radar_file_csh, start, finish, no_down=False):

# Parse radars out of the shell script - NOTE - line that contains radar list is line=6 HARDCODED

    fhandle    = open(radar_file_csh)
    all_lines  = fhandle.readlines()
    radar_list = all_lines[6].split("(")[1].split(")")[0].split()
    fhandle.close()

# Now create the calls to prep_nexrad.py
    
    for radar in radar_list:
        print(" \n Now processing %s \n " % radar)

        cmd = _prep_string % ( start.strftime("%Y,%m,%d,%H"), finish.strftime("%Y,%m,%d,%H"), "", radar, radar, "&" )
        
        if debug:
            print(cmd)

        os.system(cmd)

    old_size = 0
    new_size = get_folder_size(start_path = '.')
    cpu.sleep(_chk_dir_size)
    while old_size < new_size:
         old_size = new_size
         new_size = get_folder_size(start_path = '.')
         cpu.sleep(_chk_dir_size)

    return
#=======================================================================================================================
# getS3filelist

def getS3FileList(radar, datetime):


    noaas3 = boto3.client('s3', region_name = _region)

    files = []

    prefix = "%s/%s/" % (datetime.strftime("%Y/%m/%d"), radar)

    if debug:
        print(" \n getS3FileList string: %s \n" % prefix)

    file_list = noaas3.list_objects_v2(Bucket='noaa-nexrad-level2', Delimiter='/', Prefix=prefix)

    for i in file_list['Contents'][:]:

        if i['Key'][-2::] == 'gz' and i['Key'][-13:-11] == datetime.strftime("%H"):
            files.append(i['Key'])
        else:
            continue

    #here you can then feed the list to boto3.s3.copy_object or use wget proce

    return files

#-------------------------------------------------------------------------------
# Main program for testing...
#
if __name__ == "__main__":

#
# Command line interface 
#
    parser = OptionParser()

    parser.add_option(      "--newse",    dest="newse",    type="string", default=None, \
                                      help = "NEWSe radars description file to parse and run get_nexrad on" )

    parser.add_option("-r", "--radar",    dest="radar",    type="string", default=None, \
                                      help = "What radar to download")

    parser.add_option(      "--start",    dest="start",    type="string", default=None,  \
                                     help = "Start time of search in YYYY,MM,DD,HH")

    parser.add_option(      "--end",      dest="end",      type="string", default=None,  \
                                     help = "End time of search in YYYY,MM,DD,HH")

    parser.add_option("-d", "--dir",      dest="dir",      type="string", default=None,  \
                                     help = "directory for radar files")

    parser.add_option(      "--nthreads", dest="nthreads", type="int",    default=_nthreads, \
                                     help = "Number of download threads to run")

    (options, args) = parser.parse_args()

    if options.start == None:
        print "\n                NO START DATE AND HOUR SUPPLIED, EXITING.... \n "
        parser.print_help()
        print
        sys.exit(1)
    else:
        start   = DT.datetime.strptime(options.start, "%Y,%m,%d,%H") 

    if options.end == None:
        print "\n                NO END DATE AND HOUR SUPPLIED, EXITING.... \n "
        parser.print_help()
        print
        sys.exit(1)
    else:
        finish  = DT.datetime.strptime(options.end, "%Y,%m,%d,%H")

    if options.newse:
       print(" \n now processing NEWSe radar file....\n ")
       parse_NEWSe_radar_file(options.newse, start, finish, no_down=options.no_down)
       sys.exit(0)
        
    if options.radar == None:
        print "\n                NO RADAR SUPPLIED, EXITING.... \n "
        parser.print_help()
        print
        sys.exit(1)
    else:
        radar = options.radar
    
    if options.dir == None:
        out_dir = options.radar
    else:
        out_dir = option.dir

# Make sure we got somewhere to put the radar files
 
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

# Get file list from Amazon S3 and parse the files to get the correct hours  

    filelist = []

    ctime = start 
    
    while ctime < finish:

        newfiles = getS3FileList(radar, ctime)

        for nf in newfiles:
            filelist.append(nf)

        ctime = ctime + _delta_t

# Download each file and put it into a directory

    c0 = cpu.time()
    
    pool = Pool(processes=options.nthreads)              # set up a queue to run
    
    for file in filelist:
        
        cmd = "%s%s -P %s" % (_wget_string, file, out_dir)
            
        if debug:
            print(cmd)
            
        pool.apply_async(RunProcess, (cmd,))
    
    pool.close()
    pool.join()

    cpu0 = cpu.time() - c0
    
    print "\nDownload from Amazon took   %f  secs\n" % (round(cpu0, 3))
