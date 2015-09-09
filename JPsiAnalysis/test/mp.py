#!/usr/bin/env python
from ROOT import *
import os
#from cat_ana import *
from gen import *
import multiprocessing

def process( sample, index, file) :
  print "Making histogram %s index %d"%(sample,index)
  ana = JpsiAna( file, '%s__%s.root'%(sample,index) )
  ana.Ana()

if __name__ == "__main__" :
  name = "genJpsi"
  files = open("filelist.txt").readlines()
  p = multiprocessing.Pool(multiprocessing.cpu_count())
  path = "root://cms-xrdr.sdfarm.kr//xrd" 
  for i,file in enumerate(files) :
    p.apply_async(process, [ name, i, path+file.strip() ])
  p.close()
  p.join()
