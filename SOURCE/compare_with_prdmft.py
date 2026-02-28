#! /usr/bin/env python
import os,sys,string,math,numpy

from filecmp import dircmp
def print_diff_files(dcmp):
     for name in dcmp.diff_files:
         print("diff_file %s found in %s and %s" % (name, dcmp.left,
               dcmp.right))
     for sub_dcmp in dcmp.subdirs.values():
         print_diff_files(sub_dcmp)

dcmp = dircmp('/home/kasia/GAMMCOR/SOURCE', '/home/kasia/pr-dmft/SOURCE') 
print_diff_files(dcmp) 


#os.system("ls -l > skas")
#f=open("skas",'r')
#for line in f:
#	s=line.split(" ")
#	i=-1
#	for p in s:
#		i=i+1
#		print i,p




