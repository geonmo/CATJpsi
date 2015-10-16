#!/usr/bin/env python

import os,glob

filelist = glob.glob("TTJets_mass*.root")

for file in filelist :
  title = file.split(".root")[0]
  cmd = "./ana.py -t ttll/tree -o hist/%s %s &"%(file, file) 
  print cmd
  os.system(cmd)
