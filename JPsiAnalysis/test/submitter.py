#!/usr/bin/env python

import os,sys,glob
filelists = glob.glob("filelist*.txt")
print filelists

datasets ={}
for filelist in filelists :
  dataset = filelist.split('filelist_')[-1].split(".txt")[0]
  datasets[dataset]=filelist



for x in datasets :
  #cmd = "xrd cms-xrdr.sdfarm.kr mkdir /store/user/geonmo/catData/20151012/%s"%(x)
  #print cmd 
  cmd = "create-batch --jobName %s --maxFiles 100 --fileList %s --cfg run_TtbarDiLeptonAnalyzer.py --transferDest \"/store/user/geonmo/catData/20151012_4/%s/\""%(x,datasets[x],x)
  print cmd
  os.system(cmd)
