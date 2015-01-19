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
  def isEqual( self, jpsi1, jpsi2, exact=True ) :
    rel_pt = abs((jpsi1.pt()-jpsi2.pt())/jpsi1.pt())
    d_eta = jpsi1.eta()-jpsi2.eta()
    d_phi = TVector2.Phi_mpi_pi(jpsi1.phi()-jpsi2.phi())
    dRval = sqrt( d_eta*d_eta + d_phi*d_phi)
    if ( exact ) :
      if ( rel_pt < 1.1 and dRval < 0.025 ) :
        return True 
    else :
      #if ( rel_pt < 0.25 and rel_eta < 0.25 and rel_phi < 0.25 ) :
      if ( dRval < 0.25 ) :
        return True
    if ( exact) :
      print "Rel :  ",rel_pt, dRval
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
    nGenJpsi = TH1F("gen jpsi","n_gen_jpsi",10,0,10)
    nCatJpsi = TH1F("cat jpsi","n_cat_jpsi",10,0,10)
    nMuon = TH1F("# of muons","n_muons",10,0,10)
    trJpsi = TH1F("Trigger cat_jpsi","tr_jpsi",2,0,2)
    trJpsi.GetXaxis().SetBinLabel(1,"failed")
    trJpsi.GetXaxis().SetBinLabel(2,"passed")

    cat_vs_gen = TH2F("cat vs gen","cat_gen",10,0,10,10,0,10)
    cat_vs_gen.GetXaxis().SetTitle("# of GEN J/psi")
    cat_vs_gen.GetYaxis().SetTitle("# of CAT J/psi")

    
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
        #print "SecVtx is null."
        secVtxs=[]
        pass

      nMuon.Fill( len(muons))
      gen_jpsis =[]
      ori_gen_jpsis=[]
      """
      if ( len(muons)>0 and muons[0].pt()<4 ) :
        continue
      """
      for gen in gens :
        if ( gen.pdgId() == 443) :
          gen_jpsis.append( gen )
          ori_gen_jpsis.append(gen)
        #print gen.pdgId(),
      #print "\n" 

      cat_jpsis=[]
      for cat_jpsi in secVtxs :
        cat_jpsis.append( cat_jpsi)

      #Remove duplicate.
      for gen_jpsi1 in gen_jpsis :
        for gen_jpsi2 in gen_jpsis[ (gen_jpsis.index(gen_jpsi1)+1):] :
          if ( self.isEqual( gen_jpsi1, gen_jpsi2) ) :
            gen_jpsis.remove( gen_jpsi2)
      for cat_jpsi1 in cat_jpsis :
        for cat_jpsi2 in cat_jpsis[(cat_jpsis.index( cat_jpsi1)+1):] :
          if ( self.isEqual( cat_jpsi1, cat_jpsi2, False)) :
            cat_jpsis.remove( cat_jpsi2)

      nGenJpsi.Fill( len(gen_jpsis) )
      nCatJpsi.Fill( len(cat_jpsis) )
      cat_vs_gen.Fill( len(gen_jpsis), len(cat_jpsis))
      if ( len(cat_jpsis) > 0 ) :
        trJpsi.Fill(1)
      else :
        trJpsi.Fill(0)
      if ( len(gen_jpsis)>1 ) :
        print "Event : ",iev 
        for x in ori_gen_jpsis :
          print "GEN : ",x.pt(),x.eta(),x.phi()
        if ( len(cat_jpsis)>0 ) :
          print len(gen_jpsis), len(cat_jpsis)
          for cat_jpsi in cat_jpsis :
            print cat_jpsi.pt(), cat_jpsi.eta(), cat_jpsi.phi()


    output.Write()
    
      
   
if __name__ == "__main__" :
  filename = sys.argv[1]
  gen = Gen()
  gen.Ana(filename)
