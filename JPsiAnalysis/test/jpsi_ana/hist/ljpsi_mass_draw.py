#!/usr/bin/env python
import glob
from ROOT import *
from array import array
filelist= glob.glob("TTJets_mass*.root")

nevent = { 
            "TTJets_mass161_5_TuneZ2star_8TeV-madgraph-tauola" : 5369214,
            "TTJets_mass163_5_TuneZ2star_8TeV-madgraph-tauola" : 5365348,
            "TTJets_mass166_5_TuneZ2star_8TeV-madgraph-tauola" : 4469096,
            "TTJets_mass169_5_TuneZ2star_8TeV-madgraph-tauola" : 5202817,
            "TTJets_mass175_5_TuneZ2star_8TeV-madgraph-tauola" : 5186494,
            "TTJets_mass178_5_TuneZ2star_8TeV-madgraph-tauola" : 4733483,
            "TTJets_mass181_5_TuneZ2star_8TeV-madgraph-tauola" : 5145140,
            "TTJets_mass184_5_TuneZ2star_8TeV-madgraph-tauola" : 5249525,
        }
 
mean_gauss=[]
error_gauss = []
for file in filelist :
  title = file.split(".root")[0]
  file0 = TFile(file)
  hist1 = file0.Get("ljpsi1_mass") 
  hist2 = file0.Get("ljpsi2_mass")
 
  clone_hist = hist1.Clone()
  clone_hist2 = hist2.Clone()

  #clone_hist.Add( clone_hist2)
  clone_hist.Scale( 252.89 * 19.78 * 1000 / nevent[title])
  #clone_hist.Rebin(5)
  # Gauss
  #clone_hist.Fit("gaus")
  tf1 = TF1("f1","gaus",30,80)
  for x in range(5) :
    clone_hist.Fit(tf1)
  c1 = TCanvas("c1","c1",600,600)
  clone_hist.Draw()
  fitresult = TVirtualFitter.GetFitter()
  sig = fitresult.GetParameter(1)
  error = gitreuslt.GetErrors(1)
  mean_gauss.append(sig)
  mean_error.append(error)
  c1.SaveAs(title+"_gauss.png")

mass_various = [ 161.5, 163.5, 166.5, 169.5, 175.5, 178.5, 181.5, 184.5]
x_array = array("f", mass_various)


c1 = TCanvas("c1","c1",600,600)
y_array = array("f",mean_gauss)
h1 = TGraph( len(x_array), x_array, y_array)
h1.Draw("AL*")
c1.SaveAs("top_mass_gauss.png")
  
