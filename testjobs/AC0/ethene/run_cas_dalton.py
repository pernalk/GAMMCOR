#! /usr/bin/env python
#
import os,sys,string,math,time

DALTON_EXEC="/home/hapka/Programs/dalton/build_ifort/dalton"

# dalton file name without extension must be the argument
filename=sys.argv[1]
#
os.system(DALTON_EXEC + " -get \" AOTWOINT AOONEINT SIRIFC SIRIUS.RST rdm2.dat \" "+filename)
os.system("mv "+filename+".AOONEINT AOONEINT")
os.system("mv "+filename+".AOTWOINT AOTWOINT")
os.system("mv "+filename+".SIRIFC SIRIFC")
os.system("mv "+filename+".SIRIUS.RST SIRIUS.RST")
os.system("mv "+filename+".rdm2.dat rdm2.dat")
#
