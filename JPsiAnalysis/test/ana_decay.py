#!/usr/bin/env python
from ROOT import *
from CATJpsi.JPsiAnalysis.Lepton import Lepton
from CATJpsi.JPsiAnalysis.Muon import Muon
from CATJpsi.JPsiAnalysis.Jet import Jet
from CATJpsi.JPsiAnalysis.Electron import Electron
from CATJpsi.JPsiAnalysis.Jpsi import Jpsi

gSystem.Load("libFWCoreFWLite.so");
gSystem.Load("libDataFormatsFWLite.so");
AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,math,sys
import copy
from AbsAna import AbsAna

class DecayAna(AbsAna) :
  hist={}
  def __init__( self, infile, outfile ) :
    AbsAna.__init__(self,infile, outfile)
    self.type_list ={ 1:"GENMu", 2:"CATMu" , 3:"PFMu" , 4:"GLBMu" , 5:"TRKMu", 6:"STDMu",7:"GTrk",8:"CJpsi" }

  def BookingHist( self,output ) :
    for muon_type in self.type_list.values() : 
      self.hist[muon_type+'_muon_pt'  ] = TH1F (muon_type+"_muon_pt",   muon_type+"_muon_pt; Muon pT(GeV/c); Entries",     100,0,200)
      self.hist[muon_type+'_muon_eta' ] = TH1F (muon_type+"_muon_eta",  muon_type+"_muon_eta",                             100,-5,5)
      self.hist[muon_type+'_muon_phi' ] = TH1F (muon_type+"_muon_phi",  muon_type+"_muon_phi",                             100,-math.pi,math.pi)
      self.hist[muon_type+'_muon_mass'] = TH1F (muon_type+"_muon_mass", muon_type+"_muon_mass;Muon mass(GeV/c^{2});Entries", 100,0.08,0.12)
      self.hist[muon_type+'_muon_entries'] = TH1F (muon_type+"_nMuons", muon_type+"_nMuons", 4,0,4)
    self.hist['decayed_muon_type'] = TH1F("decayed_muon_type","Found muon's type from J/#psi;Muon type;Entries(Total)",8,1,9)
    for idx, type in enumerate(self.type_list.values()) :
      self.hist['decayed_muon_type'].GetXaxis().SetBinLabel(1+idx,type)
    self.hist['nMuons'] = TProfile("n_Muon","Number of found Muons per event",9,1,10)
    self.hist['nMuons'].GetXaxis().SetBinLabel(1,"Total")
    for idx, type in enumerate(self.type_list.values()) :
      self.hist['nMuons'].GetXaxis().SetBinLabel(idx+2,type)
  

  def FillHist( self, muon, muon_type="GENMu") :
    #print muon_type, muon.pt(), muon.eta(), muon.phi,muon.mass()
    self.hist[muon_type+"_muon_pt"].Fill( muon.pt())
    self.hist[muon_type+"_muon_eta"].Fill( muon.eta() )
    self.hist[muon_type+"_muon_phi"].Fill( muon.phi() )
    self.hist[muon_type+"_muon_mass"].Fill( muon.mass() )

    for idx, type in self.type_list.iteritems() :
      if ( muon_type == type) :
        self.hist['decayed_muon_type'].Fill( idx )
    
  def Ana( self ) :
    output = TFile(self.outfile,"RECREATE")
    self.BookingHist(outfile)

    events = self.events

    goodVTXLabel, GVTX = "goodOfflinePrimaryVertices", Handle("std::vector<reco::Vertex>")
    catMuonLabel, catMuonHandle = "catMuons", Handle("std::vector<cat::Muon>")
    genParticleLabel, genParticle = "genParticles", Handle("vector<reco::GenParticle>")
    pfAllMuonLabel, pfAllMuon = "pfAllMuonsPFlow", Handle("vector<reco::PFCandidate>")
    recoMuonLabel, recoMuon = "muons", Handle("vector<reco::Muon>")
    patTrackCandLabel, patTrackCand = "patTrackCandsPFlow", Handle("vector<pat::GenericParticle>")
    associatedpfAllMuonLabel, assopfAllMuon =  "pfAllMuonsGenParticlesMatch" , Handle("edm::Association<vector<reco::GenParticle> >")
    catJpsiLabel, catJpsiHandle = "catSecVertexs", Handle("std::vector<cat::SecVertex>")    


    for iev,event in enumerate(events):
      event.getByLabel(goodVTXLabel, GVTX)
      event.getByLabel(catMuonLabel, catMuonHandle)
      event.getByLabel(catJpsiLabel, catJpsiHandle)
      event.getByLabel(genParticleLabel, genParticle)
      event.getByLabel(pfAllMuonLabel, pfAllMuon)
      event.getByLabel(patTrackCandLabel, patTrackCand)
      event.getByLabel(associatedpfAllMuonLabel, assopfAllMuon )
      event.getByLabel(recoMuonLabel,recoMuon )
      
      
      catJpsis = None
      goodvtx = GVTX.isValid()
      try :      
        genParticles = genParticle.product()
        recoMuons = recoMuon.product()
      except RuntimeError :
        print iev 
        print "Warning! Among of one's product is null."
        continue
      try :
        pfAllMuons = pfAllMuon.product()
      except RuntimeError :
        print iev
        print "Warning! pfMuon is missing."
        pfAllMuons =[]
      try :
        catMuons = catMuonHandle.product()
      except RuntimeError :
        print "Warning! Cat Muon."
        catMuons =[]

      try : 
        patTrackCands = patTrackCand.product()
      except RuntimeError :
        print iev
        print "patTrackCand error"
        patTrackCands =[]

      try :
        catJpsis = catJpsiHandle.product()
      except RuntimeError :
        print "CAT Jpsi collection is empty or null ptr."
        catJpsis =[]

      gen_list =[]
      #if ( genParticles is None) : continue
      nMuon=[0,0,0,0,0,0,0,0,0]
      nMuon[0] +=1
      for gen in genParticles :
        if (abs(gen.pdgId())==443 and gen.numberOfDaughters()==2) :
          if ( abs(gen.daughter(0).pdgId())==13 and abs(gen.daughter(1).pdgId())==13) :
            gen_list.append( gen.daughter(0))
            gen_list.append( gen.daughter(1))
            nMuon[1] +=2
            self.FillHist( gen.daughter(0), "GENMu")
            self.FillHist( gen.daughter(1), "GENMu")
            for catJpsi in self.cleaning(catJpsis,False) :
              if self.isEqual(gen, catJpsi,False) :
                self.FillHist( gen.daughter(0), "CJpsi")
                self.FillHist( gen.daughter(1), "CJpsi")
                nMuon[8] +=2
                

      for gen in gen_list :
        for pfmuon in pfAllMuons :
          if ( self.isEqual( gen, pfmuon, False)) :
            self.FillHist( pfmuon, "PFMu")
            nMuon[3] +=1 
        for generalTrack in patTrackCands :
          if ( self.isEqual( gen, generalTrack, False)) :
            self.FillHist( generalTrack, "GTrk")
            nMuon[7] +=1 
        for catmuon in catMuons :
          if ( self.isEqual( gen, catmuon, False)) :
            self.FillHist( catmuon, "CATMu") 
            nMuon[2] +=1 
        for recomuon in recoMuons :
          if ( self.isEqual( gen, recomuon, False)) :
            if ( recomuon.isStandAloneMuon()) :
              self.FillHist( recomuon, "STDMu")
              nMuon[6] +=1 
            if ( recomuon.isTrackerMuon()) :
              self.FillHist( recomuon,"TRKMu")
              nMuon[5] +=1 
            if ( recomuon.isGlobalMuon()) :
              self.FillHist( recomuon,"GLBMu")
              nMuon[4] +=1
     
      for x in range(9)[1:] :
        self.hist['nMuons'].Fill(x+1,nMuon[x])
        self.hist[self.type_list[x]+'_muon_entries'].Fill( nMuon[x])
        


    output.Write()

  def Print(self) :
    print self.inputfile, self.outfilename

if __name__ == "__main__" :
  infile = sys.argv[1]
  outfile = sys.argv[2]
  ana = DecayAna( infile, outfile )
  ana.Ana()
