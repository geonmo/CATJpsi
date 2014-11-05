#!/usr/bin/env python
import ROOT
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle

import os

signal_path = "/cms/data/xrd/store/user/geonmo/Jpsi_20141104_MC/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola/"

class JpsiAna :
  def __init__( self, infile, mode, outfile ) :
    self.inputfile = infile
    self.mode = mode
    self.outfile = ROOT.TFile("hist/"+outfile+".root","RECREATE")
  def Ana( self ) :
    events = Events(signal_path+self.inputfile)
    goodVTXLabel, GVTX = "goodOfflinePrimaryVertices", Handle("std::vector<reco::Vertex>")
    catMuonLabel, catMuon = "catMuons", Handle("std::vector<cat::Muon>")
    for iev,event in enumerate(events): 
      event.getByLabel(goodVTXLabel, GVTX)
      event.getByLabel(catMuonLabel, catMuon)
      goodvtx = GVTX.isValid()
      print catMuon.size()

  def Print(self) :
    print self.inputfile, self.mode, self.outfile.GetName()
    
    

if __name__ == "__main__" :
  files = os.listdir(signal_path)
  for i,file in enumerate(files) :
    for mode in ["MuMu","MuEl","ElEl"] :
      ana = JpsiAna( file, mode,"hist_%s_%04d"%(mode,i) )
      ana.Ana()
