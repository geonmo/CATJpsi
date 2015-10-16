#!/usr/bin/env python

import os, sys, glob


map1={}
dirs= glob.glob("/xrootd/store/user/geonmo/catData/20151012_4/*")

for dir in dirs :
  dataset = dir.split("/")[-1]
  map1[dataset]= dir

for x in map1.keys() :
  cmd = "hadd %s.root %s/*.root"%(x, map1[x])
  print cmd
  os.system(cmd)
