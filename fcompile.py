#!/usr/bin/env python

                     #  -I/usr/local/include","-lgomp -L/usr/local/lib -lnetcdf -lnetcdff \
                     #  -L/usr/local/Cellar/zlib/1.2.8/lib -lhdf5 -lhdf5_hl -lz"], 
import sys
import os
import glob
import string
from optparse import OptionParser

f2py_only    = ["cressman_kdtree.f90"]

fortran_only = ["kdtree2.f90"]

preprocess   = ""
preprocess   = "-DF2PY_REPORT_ON_ARRAY_COPY"

fopts        = {'pgf': ['pg',"-Bdynamic -fast -mp",""], \
                'gnu': ["gnu95","-O3 -funroll-loops -fopenmp -fPIC -I/usr/local/include", "-lgomp -L/usr/local/lib" ],
                'intel': ['intelem',"-v -O3", "-liomp5 -L/opt/local/netcdf4m64/lib -lnetcdf -lnetcdff \
                          -L/opt/local/hdf5m64/lib -lhdf5 -lhdf5_hl -lz"]}

parser = OptionParser()
parser.add_option("-f","--fc",dest="compiler",type="string", default='gnu', \
                  help = "fortran compiler to be used valid compiler: [gfortran, intel, ppgf, default is gfortran]")

(options, args) = parser.parse_args()

objects = ""

# remove all module files in directory - this can trip you up bad!

cmd = 'rm *.mod'
print "\n=====================================================\n"
print("   ---> Removing all module files...safety first!")
print "\n=====================================================\n"
ret = os.system(cmd)

# First compiler the fortran only files
 
for item in fortran_only:
    print "\n=====================================================\n"
    print "  Compiling file: %s using the %s compiler" % (item, options.compiler)
    print "\n====================================================="
    if options.compiler == 'gnu':
        print("\n  Using GNU gfortran compiler \n")
        cmd = 'gfortran %s -c %s ' % (fopts['gnu'][1],item)
    if options.compiler == 'intel':
        print("\n  Using Intel compiler \n")
        cmd = 'ifort %s -c %s ' % (fopts['intel'][1],item)
    if options.compiler == 'pgf':
        print("\n  Using Portland group compiler \n")
        cmd = 'pgf90 %s -c %s ' % (fopts['pgf'][1],item)

    print("  "+cmd+"\n")
    ret = os.system(cmd)

    if ret == 0:
        objects = objects + " %s.o" % item.split(".")[0]

print "\nFortran object files compiled:  %s \n" % objects    
print "\n=====================================================\n"

for item in f2py_only:
    ret = os.system('rm %s.so ' % (item.split(".")[0]))

    if options.compiler == 'gnu':

        cmd = ("f2py --debug --fcompiler=%s --f90flags='%s' %s -c -m %s %s %s %s " % (fopts['gnu'][0],fopts['gnu'][1], \
                preprocess, item.split(".")[0], item, objects, fopts['gnu'][2]))

    if options.compiler == 'intel':

        cmd = ('f2py --fcompiler="intelem" "--f90flags= -nofor-main %s" %s -c -m %s %s %s %s' % (fopts['intel'][0], \
               fopts['intel'][1], preprocess, item.split(".")[0], item, objects,fopts['intel'][2]))

    if options.compiler == 'pgf':
        cmd = ('f2py --fcompiler=%s "--f90flags=%s" %s -c -m %s %s %s %s' % (fopts['pgf'][0],fopts['pgf'][1], \
               preprocess, item.split(".")[0], item, objects,fopts['pgf'][2]))

    print("\n\n ==> %s \n\n" % cmd)

    ret = os.system(cmd)

    print "\n==========================================================================================\n"

    if ret == 0:
        print("\nSuccessfully compiled file: %s " % item)
        print("\nObject file is: %s " % (item.split(".")[0] + ".so"))
    else:
        print "   ERROR !!!!!   ERROR--> unsuccessful compile file: %s\n" % item
        sys.exit(0)

    cmd = "python -c 'import %s'" % item.split(".")[0]
    ret = os.system(cmd)

    if ret == 0:
        print("\n   --> %s was successfully imported into python, you should be good to go...\n" % item)
    else:
        print("\n\n!!!!!! fcompile ERROR, compile must have failed as %s cannot be loaded...." % item)
        print("\n error: %s" % ret)

    print("\n==========================================================================================\n")
