#!/usr/bin/env python
from ROOT import *
import os
from cat_ana import *
import multiprocessing

def process( sample, index, file) :
  print "Making histogram %s index %d"%(sample,index)
  ana = JpsiAna( file, '%s__%s'%(sample,index) )
  ana.Ana()

if __name__ == "__main__" :
  name = "FullLepton"
  files = open("filelist.txt").readlines()
  p = multiprocessing.Pool(multiprocessing.cpu_count())
  path = "/cms/data/xrd" 
  for i,file in enumerate(files) :
    p.apply_async(process, [ name, i, path+file.strip() ])
  p.close()
  p.join()
