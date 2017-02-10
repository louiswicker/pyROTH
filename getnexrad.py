import datetime as DT
import sys, os
import boto3
from optparse import OptionParser
import datetime as DT
from multiprocessing import Pool
import time as cpu


_region        = 'us-east-1'
_delta_t       = DT.timedelta(hours=1)  
_wget_string   = "wget https://noaa-nexrad-level2.s3.amazonaws.com/"
_pyRoth_string = "pyROTH.py -d %s -u None -w -o roth_%s >& log_%s"
_pyDart_string = 'pyDart.py -d roth_%s "*VR.out" --ascii2hdf >& log_%s'
_merge_string  = 'pyDart.py -f %s -d roth_%s "*VR.h5" --merge >& log_%s'

debug = True
_nthreads = 2

#=======================================================================================================================
# RunMember is a function that runs a system command

def RunMember(cmd):
    print("\n Executing command:  %s " % cmd)
    os.system(cmd)
    print("\n %s completed...." % cmd)
    return


def getfilelist(radar, datetime):


    noaas3 = boto3.client('s3', region_name = _region)

    files = []

    prefix = "%s/%s/" % (datetime.strftime("%Y/%m/%d"), radar)

    file_list = noaas3.list_objects_v2(Bucket='noaa-nexrad-level2', Delimiter='/', Prefix=prefix)

    for i in file_list['Contents'][:]:

#HERE YOU CAN SUBSET THE FULL LIST of radar key which looks like: HERE I just grab the 23utc files.  

# '2016/05/23/KTLX/KTLX20110427_000157_V03.gz'

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
    parser.add_option("-r", "--radar",       dest="radar",     default=None,  type="string", help = "What radar to download")
    parser.add_option(      "--start",       dest="start",     type="string", help = "Start time of search in YYYY,MM,DD,HH")
    parser.add_option(      "--end",         dest="end",       type="string", help = "End time of search in YYYY,MM,DD,HH,MM")
    parser.add_option("-d", "--dir",         dest="dir",       default=None, type="string", help = "directory for radar files")
    parser.add_option(      "--nthreads", dest="nthreads", type="int",  default=_nthreads, help = "Number of threads to run")

    (options, args) = parser.parse_args()

    start   = DT.datetime.strptime(options.start, "%Y,%m,%d,%H") 

    finish  = DT.datetime.strptime(options.end, "%Y,%m,%d,%H")

    radar = options.radar
    
    if options.dir == None:
        out_dir = options.radar
    else:
        out_dir = option.dir

 
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

# Get file list from Amazon S3   
    filelist = []
    
    while start < finish:
        begin = start 

        newfiles = getfilelist(radar, begin)
        for nf in newfiles:
            filelist.append(nf)

        start = start + _delta_t

# Download each file and put it into a directory
# 
#     c0 = cpu.time()
#     
#     pool = Pool(processes=options.nthreads)              # set up a queue to run
#     
#     for file in filelist:
#     
#         cmd0 = "%s%s" % (_wget_string, file)
#             
#         cmd1 = "mv %s %s/%s" % (file[-26:], out_dir, file[-26:])
# 
#         cmd = "%s ; %s " % (cmd0, cmd1)
#         
#         if debug:
#             print(cmd)
#             
#         pool.apply_async(RunMember, (cmd,))
#     
# 
#     pool.close()
#     pool.join()
# 
#     cpu0 = cpu.time() - c0
#     
#     print "\nDownload from Amazon took   %f  secs\n" % (round(cpu0, 3))

# Process the radar files threaded

    c0 = cpu.time()
    
    cmd = _pyRoth_string % ( radar, radar, radar )
    if debug:
        print(cmd)
    os.system(cmd)

    cmd = _pyDart_string % ( radar, radar )
    if debug:
        print(cmd)
    os.system(cmd)
    
    VR_file = "obs_seq_%s_%s_VR.h5" % (start.strftime("%Y_%m_%d"), radar)
    cmd = _merge_string % ( VR_file, radar, radar )
    if debug:
        print(cmd)
        
    os.system(cmd)
    
    cpu0 = cpu.time() - c0
    
    print "\nProcessing %s files took   %f  secs\n" % (radar, round(cpu0, 3))
   
