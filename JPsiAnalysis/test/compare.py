#!/usr/bin/env python
from ROOT import *
gSystem.Load("libFWCoreFWLite.so");
gSystem.Load("libDataFormatsFWLite.so");
AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,math,sys




def hist_maker(name, title, bin_set, x_name, y_name, tr, br, cut,color=kBlack):
  hist = TH1F(name, title, bin_set[0], bin_set[1], bin_set[2])
  hist.GetXaxis().SetTitle(x_name)
  hist.GetYaxis().SetTitle(y_name)
  hist.SetLineColor(color)
  #hist.Sumw2()
  hist.SetStats(0)
  tr.Project(name, br, cut)
  return hist

def bin_mod(br, name_tag, tr, cut):
  tr.Draw(br+">>h_"+name_tag, cut)
  tmp = gDirectory.Get("h_"+name_tag)
  bin_max = tmp.GetXaxis().GetXmax()
  bin_min = tmp.GetXaxis().GetXmin()
  bin_d = bin_max - bin_min
  if bin_d >200:
    return [200, bin_min, bin_max]
  else:
    return [100, bin_min, bin_max]


file = TFile.Open(sys.argv[1])


matched = file.Get("matched_jpsis")
unmatched = file.Get("unmatched_jpsis")

variable = ["pt","eta","phi","mass","lxy","l3D","vProb","sigmalxy","dca","cxPtHypot","cxPtAbs"]

"""
bins = bin_mod( "pt","hhhhh",matched, "1")
h1 = hist_maker( variable[0],variable[0], bins, "X","Y",matched, variable[0],"1",kRed)
h1.Draw()
"""
for i,x in enumerate(variable) :
  bins = bin_mod(x,x+"temp",matched,"1")
  canvas_name = "cc%d"%(i) 
  c = TCanvas(canvas_name,canvas_name,800,600)
  h1 = hist_maker( variable[i], variable[i]+"_matched",bins, variable[i],"Y",matched, variable[i], "1",kBlue)
  h2 = hist_maker( variable[i]+"_1", variable[i]+"_unmatched",bins, variable[i],"Y",unmatched, variable[i], "1",kRed)
  h1.Scale( 1/h1.GetEntries())
  h2.Scale( 1/h2.GetEntries())
  h1.Draw()
  h2.Draw("same")
  leg = c.BuildLegend()
  h1.SetTitle(variable[i])
  c.SaveAs(canvas_name+".png")
