#!/usr/bin/env python
from ROOT import *
import glob

info = { "TTJets_FullLeptMGDecays_8TeV-madgraph-tauola": (26.54, 12011428) , "TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola": (110.778, 31004465), "TTJets_HadronicMGDecays_8TeV-madgraph":(115.564,10537444), "DYJetsToLL_M-10To50filter_8TeV-madgraph":(907.3155,7132223),  "DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball":(3503.71,30459503), "T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola":(11.1773,497658),"Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola":(11.1773,493460), "WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball":(37509,57709905), "WW_TuneZ2star_8TeV_pythia6_tauola":(56,10000431), "WZ_TuneZ2star_8TeV_pythia6_tauola":(33.6,10000283), "ZZ_TuneZ2star_8TeV_pythia6_tauola":(7.6,9799908)}

RD_label = ["DoubleMuParked","DoubleElectron","MuEG"]
filelist = glob.glob("*.root")
mc_list=[]
rd_list=[]


for file in filelist :
  if( file[:-5] in RD_label ) : rd_list.append(file)
  elif ( file.find("mass") != -1) : continue
  else : mc_list.append(file)

cut_flow =[]
num_event = {}


for rd in rd_list:
  file0 = TFile(rd)
  c1 = TCanvas(rd)
  tdir = file0.Get("ttll/tree")
  title = rd.split(".root")[0]
  num_event[title] ={}
  for i in range(3) :
    tdir.Draw("jpsi_mass >>hnew","step==5 && channel == %d && jpsi_mva > -0.0360"%(i) )
    hnew = file0.Get("hnew")
    num_event[title][i] = (hnew.GetEntries())
  c1.SaveAs(title+"_num_ofJpsi.png")

for key in num_event.keys() :
  print key, 
  for i in range(3) :
    print " , ",num_event[key][i],
  print " ,"

    
print "total , ",
for i in range(3) :
  sum=0
  for key in num_event.keys() :
    sum = sum+ num_event[key][i]
  print sum," , ",
print

