#!/usr/bin/env python
from ROOT import *
gSystem.Load("libFWCoreFWLite.so");
gSystem.Load("libDataFormatsFWLite.so");
AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,math,sys

class Gen :
  def __init__(self) :
    pass
  def Ana(self,infile) :
    file = TFile(infile)
    print "Load file : "+file.GetName()
    if ( file.IsZombie()) :
      print "%s is corruct!"%(infile)
      return
    file.Close()

    events = Events(infile)
    goodVTXLabel, GVTX = "goodOfflinePrimaryVertices", Handle("std::vector<reco::Vertex>")
    catMuonLabel, catMuon = "catMuons", Handle("std::vector<cat::Muon>")
    catElectronLabel, catElectron = "catElectrons", Handle("std::vector<cat::Electron>")
    catSecVertexLabel, catSecVtx = "catSecVertexs", Handle("std::vector<cat::SecVertex>")
    catGenLabel, genParticle = "genParticles", Handle("std::vector<reco::GenParticle>")
    catJetLabel, catJet = "catJets", Handle("std::vector<cat::Jet>")
    
    for iev,event in enumerate(events):
      #print iev 
      event.getByLabel(goodVTXLabel, GVTX)
      event.getByLabel(catMuonLabel, catMuon)
      event.getByLabel(catElectronLabel, catElectron)
      event.getByLabel(catSecVertexLabel, catSecVtx)
      event.getByLabel(catJetLabel, catJet)
      event.getByLabel(catGenLabel, genParticle)

      goodvtx = GVTX.isValid()
      secVtxs = []
      try :      
        muons = catMuon.product()
        electrons = catElectron.product()
        jets = catJet.product()
      except RuntimeError :
        print "Warning! Among of one's product is null."
        continue

      try :      
        secVtxs = catSecVtx.product()
      except RuntimeError :
        print "SecVtx is null."
        secVtxs=[]
        pass
   
if __name__ == "__main__" :
  filename = sys.argv[1]
  gen = Gen()
  gen.Ana(filename)
