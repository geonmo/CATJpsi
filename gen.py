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
      if ( rel_pt < 1e-3 and dRval < 1e-3 ) :
        return True 
    else :
      if ( rel_pt < 0.05 and dRval < 0.05 ) :
        return True
    return False
   
  def cleaning( self, jpsis,exact=True) :
    c_jpsis =[]
    for jpsi in jpsis :
      duplicate = False
      for c_jpsi in c_jpsis :
        if ( self.isEqual (jpsi,c_jpsi,exact)) :
          duplicate = True
      if ( not duplicate ) :
        c_jpsis.append(jpsi)
    return c_jpsis

  def FillNtuple(self, jpsi_tuple, jpsiList) :
    for jpsi in jpsiList :
      jpsi_tuple.Fill( jpsi.pt(), jpsi.eta(), jpsi.phi(),jpsi.mass(),0.0,0.0,jpsi.lxy(),jpsi.l3D(),jpsi.vProb(),jpsi.sigmalxy(),jpsi.dca(),jpsi.cxPtHypot(),jpsi.cxPtAbs());
  """
  def SaveNtuple(self, filename, jpsi_tuple) :
    file = TFile(filename, "RECREATE")
    jpsi_tuple.Write()
    file.Close()
  """

     
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


    #output = TFile("outputJpsi.root","RECREATE")
    #nGenJpsi = TH1F("n_gen_jpsi","gen jpsi",10,0,10)
    #nCatJpsi = TH1F("n_cat_jpsi","cat jpsi",10,0,10)
    #nMuon = TH1F("n_muons","# of muons",10,0,10)
    #trJpsi = TH1F("tr_jpsi","Trigger cat_jpsi",2,0,2)
    #trJpsi.GetXaxis().SetBinLabel(1,"failed")
    #trJpsi.GetXaxis().SetBinLabel(2,"passed")

    #cat_vs_gen = TH2F("cat_gen","cat vs gen",10,0,10,10,0,10)
    #cat_vs_gen.GetXaxis().SetTitle("# of GEN J/psi")
    #cat_vs_gen.GetYaxis().SetTitle("# of CAT J/psi")
    

    output = TFile("output.root","RECREATE")

    #matched_jpsi_tuple = TNtuple("matched_jpsis","Matched J/Psi's information","pt:eta:phi:mass:chi2:ndof:lxy:l3D:vProb:sigmalxy:dca:cxPtHypot:cxPtAbs")
    matched_jpsis_ntuple = TNtuple("matched_jpsis","Matched J/Psi's information","pt:eta:phi:mass:chi2:ndof:lxy:l3D:vProb:sigmalxy:dca:cxPtHypot:cxPtAbs")
    unmatched_jpsis_ntuple = TNtuple("unmatched_jpsis","UnMatched J/Psi's information","pt:eta:phi:mass:chi2:ndof:lxy:l3D:vProb:sigmalxy:dca:cxPtHypot:cxPtAbs")
    #matched_jpsis_ntuple = TNtuple("matched_jpsis","Matched J/Psi's information","pt:eta:phi:mass")
    #unmatched_jpsis_ntuple = TNtuple("unmatched_jpsis","Un-Matched J/Psi's information","pt:eta:phi:mass")
    found_jpsi  = TH1F("isMatchedJpsi","Matched J/#psi between gen and cat.",2,0,2)
    found_jpsi.GetXaxis().SetBinLabel(1,"Failed")
    found_jpsi.GetXaxis().SetBinLabel(2,"Successed")
    cleaned_found_jpsi  = TH1F("isMatchedJpsi_cleaned","Cleaned matched J/#psi between gen and cat if only J/#psi to mu+mu-.",2,0,2)
    cleaned_found_jpsi.GetXaxis().SetBinLabel(1,"Failed")
    cleaned_found_jpsi.GetXaxis().SetBinLabel(2,"Successed")
    nJpsi = TNtuple("nJpsi","Number of J/#psi s","nJpsi")
    
    for iev,event in enumerate(events):
      #if iev > 300 : break
      #event.getByLabel(goodVTXLabel, GVTX)
      event.getByLabel(catMuonLabel, catMuon)
      #event.getByLabel(catElectronLabel, catElectron)
      event.getByLabel(catSecVertexLabel, catSecVtx)
      #event.getByLabel(catJetLabel, catJet)
      event.getByLabel(catGenLabel, genParticle)

      goodvtx = GVTX.isValid()
      secVtxs = []
      try :      
        muons = catMuon.product()
        #electrons = catElectron.product()
        #jets = catJet.product()
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

      #nMuon.Fill( len(muons))
      gen_jpsis =[]
      ori_gen_jpsis=[]
      c_gen_jpsis=[]
      for gen in gens :
        if ( gen.pdgId() == 443 and gen.pt()>1.0 and abs(gen.eta())<2.5  and gen.mass()>3 and gen.mass()<3.2) :
          if ( gen.numberOfDaughters() ==2 ) :
            if ( abs(gen.daughter(0).pdgId())==13 and abs(gen.daughter(1).pdgId())==13 ) :
              gen_jpsis.append( gen )
              #print "GEN Found"
      c_gen_jpsis = self.cleaning( gen_jpsis,True)
      

      cat_jpsis=[]
      for cat_jpsi in secVtxs :
        if ( cat_jpsi.pt()>1.0 and abs(cat_jpsi.eta())<2.5 and cat_jpsi.mass()>3 and cat_jpsi.mass()<3.2 ) :
          cat_jpsis.append( cat_jpsi)
      c_reco_jpsis = self.cleaning( cat_jpsis,False)
      
      if ( len( c_gen_jpsis) >0 ) :
        print "evt : %d //  # of gen J/psi  : %d, # of reco J/psi : %d\t"%(iev, len(c_gen_jpsis), len(c_reco_jpsis))
      #Remove duplicate.
      matched_jpsis =[]
      unmatched_jpsis =[]
      
      for reco_jpsi in c_reco_jpsis :
        matching= False
        for gen_jpsi in c_gen_jpsis :
          if ( self.isEqual( reco_jpsi, gen_jpsi,False) ) :
            matching = True
        if ( matching ) :
            matched_jpsis.append( reco_jpsi)
        else :
            unmatched_jpsis.append( reco_jpsi)

      if ( len(matched_jpsis)>0 ) :
        print "At least one jpsis is matched. ( %d)"%(len(matched_jpsis))
        found_jpsi.Fill(1)
        if ( len(c_gen_jpsis) >0 ) : cleaned_found_jpsi.Fill(1)
        nJpsi.Fill( len(c_reco_jpsis))
      else :
        #print "No jpsi is matched."
        if ( len(c_gen_jpsis) >0 ) : cleaned_found_jpsi.Fill(0)
        found_jpsi.Fill(0)
        nJpsi.Fill( 0 )

      self.FillNtuple(matched_jpsis_ntuple,matched_jpsis)
      self.FillNtuple(unmatched_jpsis_ntuple,unmatched_jpsis)

    output.Write()
    output.Close()
    
      
   
if __name__ == "__main__" :
  filename = sys.argv[1]
  gen = Gen()
  gen.Ana(filename)
