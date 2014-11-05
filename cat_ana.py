#!/usr/bin/env python
import ROOT
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle

import os

class JpsiAna :
  def __init__( self, infile, mode, outfile ) :
    self.inputfile = infile
    self.mode = mode
    self.outfile = ROOT.TFile("hist/"+outfile+".root","RECREATE")
  def Ana( self ) :
    pass
  def Print(self) :
    print self.inputfile, self.mode, self.outfile.GetName()
    
    

if __name__ == "__main__" :
  signal_path = "/cms/data/xrd/store/user/geonmo/Jpsi_20141104_MC/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola/"
  files = os.listdir(signal_path)
  for i,file in enumerate(files) :
    for mode in ["MuMu","MuEl","ElEl"] :
      ana = JpsiAna( file, mode,"hist_%s_%04d"%(mode,i) )
      ana.Print()
