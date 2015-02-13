#!/usr/bin/env python
from ROOT import *
gSystem.Load("libFWCoreFWLite.so");
gSystem.Load("libDataFormatsFWLite.so");
AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,math,sys

class MotherTracking :
  def __init__(self) :
    pass
  def isFromB( self, particle ) :
    print "pdgId : ",particle.pdgId()
    if (particle is None ) :
      return False
    if ( abs(particle.pdgId()) == 5 ) :
      return True
    flag = False
    nMother = particle.numberOfMothers()
    for idx in range(nMother) :
      if ( self.isFromB( particle.mother(idx) ) ) :
        flag = True
    return flag

  #def isFromTop( self, particle ) :
    

  def isEqual( self, jpsi1, jpsi2, exact=True ) :
    rel_pt = abs((jpsi1.pt()-jpsi2.pt())/jpsi1.pt())
    d_eta = jpsi1.eta()-jpsi2.eta()
    d_phi = TVector2.Phi_mpi_pi(jpsi1.phi()-jpsi2.phi())
    dRval = sqrt( d_eta*d_eta + d_phi*d_phi)
    if ( exact ) :
      if ( rel_pt < 1 and dRval < 0.025 ) :
        return True 
    else :
      if ( rel_pt < 0.25 and rel_eta < 0.25 and rel_phi < 0.25 ) :
        return True
    return False
 
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


    output = TFile("outputJpsi.root","RECREATE")

    for iev,event in enumerate(events):
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
        gens = genParticle.product()
      except RuntimeError :
        print "Warning! Among of one's product is null."
        continue

      try :      
        secVtxs = catSecVtx.product()
      except RuntimeError :
        secVtxs=[]
        pass

      gen_jpsis =[]
      for gen in gens :
        if ( gen.pdgId() == 443 ) :
          gen_jpsis.append( gen )

      for jpsi in gen_jpsis :
        if ( self.isFromB( jpsi ) ) :
          print "it was come from B"
        else :
          print "it did not be come from B"
          


    output.Write()
    
      
   
if __name__ == "__main__" :
  filename = sys.argv[1]
  gen = MotherTracking()
  gen.Ana(filename)
