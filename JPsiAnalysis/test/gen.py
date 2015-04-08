#!/usr/bin/env python
from ROOT import *
gSystem.Load("libFWCoreFWLite.so");
gSystem.Load("libDataFormatsFWLite.so");
AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,math,sys
from TopAna import TopAna 

class JpsiAna(TopAna) : 
  def __init__(self,infile, outfile) :
    TopAna.__init__(self,infile,outfile)

  def FillNtuple(self, jpsi_tuple, jpsiList, gen_jpsis) :
    for jpsi,match_list in jpsiList :
      if len(match_list) ==0 :
        jpsi_tuple.Fill( jpsi.pt(), jpsi.eta(), jpsi.phi(),jpsi.mass(),0.0,0.0,jpsi.lxy(),jpsi.l3D(),jpsi.vProb(),jpsi.sigmalxy(),jpsi.dca(),jpsi.cxPtHypot(),jpsi.cxPtAbs(),0,0);
      else :
        target_gen = gen_jpsis[match_list[0]]
        isFromB = self.isFromB( target_gen )
        isFromT = self.isFromTop( target_gen) or self.isFromTopBar( target_gen) 
        jpsi_tuple.Fill( jpsi.pt(), jpsi.eta(), jpsi.phi(),jpsi.mass(),0.0,0.0,jpsi.lxy(),jpsi.l3D(),jpsi.vProb(),jpsi.sigmalxy(),jpsi.dca(),jpsi.cxPtHypot(),jpsi.cxPtAbs(),isFromB,isFromT)

  def Ana(self) :
    events = self.events
    goodVTXLabel, GVTX = "goodOfflinePrimaryVertices", Handle("std::vector<reco::Vertex>")
    catMuonLabel, catMuon = "catMuons", Handle("std::vector<cat::Muon>")
    catElectronLabel, catElectron = "catElectrons", Handle("std::vector<cat::Electron>")
    catSecVertexLabel, catSecVtx = "catSecVertexs", Handle("std::vector<cat::SecVertex>")
    catGenLabel, genParticle = "genParticles", Handle("std::vector<reco::GenParticle>")
    catJetLabel, catJet = "catJets", Handle("std::vector<cat::Jet>")


    output = TFile(self.outfile,"RECREATE")

    matched_jpsis_ntuple = TNtuple("matched_jpsis","Matched J/Psi's information","pt:eta:phi:mass:chi2:ndof:lxy:l3D:vProb:sigmalxy:dca:cxPtHypot:cxPtAbs:isFromB:isFromTop")
    unmatched_jpsis_ntuple = TNtuple("unmatched_jpsis","UnMatched J/Psi's information","pt:eta:phi:mass:chi2:ndof:lxy:l3D:vProb:sigmalxy:dca:cxPtHypot:cxPtAbs:isFromB:isFromTop")
    found_jpsi  = TH1F("isMatchedJpsi","Matched J/#psi between gen and cat.",2,0,2)
    found_jpsi.GetXaxis().SetBinLabel(1,"Failed")
    found_jpsi.GetXaxis().SetBinLabel(2,"Successed")
    cleaned_found_jpsi  = TH1F("isMatchedJpsi_cleaned","Cleaned matched J/#psi between gen and cat if only J/#psi to mu+mu-.",2,0,2)
    cleaned_found_jpsi.GetXaxis().SetBinLabel(1,"Failed")
    cleaned_found_jpsi.GetXaxis().SetBinLabel(2,"Successed")
    nGENJpsi = TNtuple("nGENJpsi","Number of GEN J/#psi s","nGENJpsi")
    nCATJpsi = TNtuple("nCATJpsi","Number of CAT J/#psi s","nCATJpsi")
    nMuon = TNtuple("nMuon","Number of Muons","nMuon")
    nMuonPair = TNtuple("nMuonPair","Number of Muon Pairs","nMuonPair")

    isRD = False


    jpsi_minPt = 4.0

    
    for iev,event in enumerate(events):
      event.getByLabel(catMuonLabel, catMuon)
      event.getByLabel(catSecVertexLabel, catSecVtx)

      secVtxs = []
      try :      
        muons = catMuon.product()
        #electrons = catElectron.product()
      except RuntimeError :
        print "Warning! CAT Muon product is null."
        continue
      try :
        event.getByLabel(catGenLabel, genParticle)
        gens = genParticle.product()
      except (ValueError, RuntimeError):
        if ( not isRD ) :
          print "No genparticle. We assume these are RD"
        isRD = True
        gens = []

      try :      
        secVtxs = catSecVtx.product()
      except RuntimeError :
        secVtxs=[]
        pass

      gen_jpsis =[]
      ori_gen_jpsis=[]
      c_gen_jpsis=[]
      for gen in gens :
        #if ( gen.pdgId() == 443 ) :
        if ( gen.pdgId() == 443 and gen.pt()>1.0 and abs(gen.eta())<2.5  and gen.mass()>3 and gen.mass()<3.2) :
          if ( gen.numberOfDaughters() ==2 ) :
            if ( abs(gen.daughter(0).pdgId())==13 and abs(gen.daughter(1).pdgId())==13 ) :
              gen_jpsis.append( gen )
      c_gen_jpsis = self.cleaning( gen_jpsis,True)
      nGENJpsi.Fill( len(c_gen_jpsis))
      

      cat_jpsis=[]
      for cat_jpsi in secVtxs :
        #if (True) :
        #if ( cat_jpsi.pt()>1.0 and abs(cat_jpsi.eta())<2.5 and cat_jpsi.mass()>3 and cat_jpsi.mass()<3.2 and cat_jpsi.vProb()>0.01 ) :
        if ( cat_jpsi.pt()>jpsi_minPt and abs(cat_jpsi.eta())<2.5 and cat_jpsi.mass()>3 and cat_jpsi.mass()<3.2 and cat_jpsi.vProb()>0.01 and cat_jpsi.dca()<0.05 and cat_jpsi.l3D()>0.2 ) :
          cat_jpsis.append( cat_jpsi)
      c_reco_jpsis = self.cleaning( cat_jpsis,False)
      nCATJpsi.Fill( len(c_reco_jpsis))
      
      matched_jpsis =[]
      unmatched_jpsis =[]

      c_muons = []
      for muon in muons :
        if ( muon.pt()>jpsi_minPt and abs(muon.eta())<2.5) :
          c_muons.append( muon)
      nMuon.Fill( len( c_muons))
      pair = 0
      for idx,muon1 in enumerate(c_muons) :
        for muon2 in c_muons[idx:] :
          if ( muon1.charge()*muon2.charge() == -1 ) :
            pair = pair +1
      nMuonPair.Fill(pair)


      if ( len( c_gen_jpsis) >0 ) :
        print "evt : %d //  # of gen J/psi  : %d, # of reco J/psi : (%d/%d), # of Muons : %d \n"%(iev, len(c_gen_jpsis), len(c_reco_jpsis), len(secVtxs), len(muons))
        for gen in c_gen_jpsis :
          print "GEN  pT : %f, Eta: %f, Phi: %f, Mass: %f\n"%(gen.pt(),gen.eta(),gen.phi(), gen.mass())
        for cat_jpsi in c_reco_jpsis :
          print "RECO pT : %f, Eta: %f, Phi: %f, Mass: %f, vProb: %f, DCA: %f, muonID %d, trackQuality : %d\n"%(cat_jpsi.pt(),cat_jpsi.eta(),cat_jpsi.phi(), cat_jpsi.mass(),cat_jpsi.vProb(), cat_jpsi.dca(),cat_jpsi.muID(), cat_jpsi.trackQuality())

      for reco_jpsi in c_reco_jpsis :
        matching= False
        match_pair = [] 
        for idx,gen_jpsi in enumerate(c_gen_jpsis) :
          if ( self.isEqual( reco_jpsi, gen_jpsi,False) ) :
            matching = True
            match_pair.append(idx)
        if ( matching ) :
            matched_jpsis.append( [reco_jpsi,match_pair])
        else :
            if ( not isRD and c_gen_jpsis>0 ) :
              unmatched_jpsis.append( [reco_jpsi,match_pair])
            elif (isRD) :
              unmatched_jpsis.append( [reco_jpsi,match_pair])
            else :  ## no genJ/psi on MC  
              continue

      if ( len(matched_jpsis)>0 ) :
        print "At least one jpsis is matched. ( %d)"%(len(matched_jpsis))
        found_jpsi.Fill(1)
        if ( len(c_gen_jpsis) >0 ) : cleaned_found_jpsi.Fill(1)
      else :
        #print "No jpsi is matched."
        if ( len(c_gen_jpsis) >0 ) : cleaned_found_jpsi.Fill(0)
        found_jpsi.Fill(0)

      self.FillNtuple(matched_jpsis_ntuple,matched_jpsis, c_gen_jpsis)
      self.FillNtuple(unmatched_jpsis_ntuple,unmatched_jpsis, c_gen_jpsis)

    output.Write()
    output.Close()
    
      
   
if __name__ == "__main__" :
  infile = sys.argv[1]
  outfile = sys.argv[2]
  gen = JpsiAna(infile, outfile)
  gen.Ana()
