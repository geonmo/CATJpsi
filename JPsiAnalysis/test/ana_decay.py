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

  def BookingHist( self,output ) :
    for muon_type in ['gen','cat','pf','global','tracker'] :
      self.hist[muon_type+'_muon_pt'  ] = TH1F (muon_type+"_muon_pt",   muon_type+"_muon_pt; Muon pT(GeV/c); Entries",     100,0,200)
      self.hist[muon_type+'_muon_eta' ] = TH1F (muon_type+"_muon_eta",  muon_type+"_muon_eta",                             100,-5,5)
      self.hist[muon_type+'_muon_phi' ] = TH1F (muon_type+"_muon_phi",  muon_type+"_muon_phi",                             100,-math.pi,math.pi)
      self.hist[muon_type+'_muon_mass'] = TH1F (muon_type+"_muon_mass", muon_type+"_muon_mass;Muon mass(GeV/c^2);Entries", 100,0.1,0.11)

  def FillHist( self, muon, muon_type="gen") :
    self.hist[(muon_type+"_muon_pt")].Fill( muon.Pt())
    self.hist[(muon_type+"_muon_eta")].Fill( muon.Eta() )
    self.hist[(muon_type+"_muon_phi")].Fill( muon.Phi() )
    self.hist[(muon_type+"_muon_mass")].Fill( muon.M() )
    
  def Ana( self ) :
    file = TFile.Open(self.inputfile)
    print "Load file : "+file.GetName()
    if ( file.IsZombie()) :
      print "%s is corruct!"%(self.inputfile)
      return
    file.Close()

    outfile = TFile(self.outfilename,"RECREATE")
    for mode in c_modes :
      for x in steps :
        outfile.mkdir("%s_s%d"%(mode,x))
    self.BookingHist(outfile)

    events = Events(self.inputfile)

    goodVTXLabel, GVTX = "goodOfflinePrimaryVertices", Handle("std::vector<reco::Vertex>")
    catMuonLabel, catMuon = "catMuons", Handle("std::vector<cat::Muon>")
    for iev,event in enumerate(events):
      event.getByLabel(goodVTXLabel, GVTX)
      event.getByLabel(catMuonLabel, catMuon)
      event.getByLabel(catElectronLabel, catElectron)
      event.getByLabel(catSecVertexLabel, catSecVtx)
      event.getByLabel(catJetLabel, catJet)

      goodvtx = GVTX.isValid()
      secVtx = None
      try :      
        muons = catMuon.product()
        electrons = catElectron.product()
        jets = catJet.product()
      except RuntimeError :
        print iev 
        print "Warning! Among of one's product is null."
        continue
      if ( len(muons)<2 and len(electrons) <2) :
        continue
      try :
        secVtxs = catSecVtx.product()
      except RuntimeError :
        print "SecVtx is missing"
        secVtxs =[]
      #print "muon : %d electorns : %d jets : %d"%(len(muons),len(electrons),len(jets))

      muon_list =[]
      elec_list =[]
      jet_list  =[]
      jpsi_list =[]
      l1jpsi_list=[]
      l2jpsi_list=[]

      #Acquiring POGs.
      for muon in muons :
        t_muon = Muon(muon,'m')
        if ( t_muon.valid() ) :
          muon_list.append( t_muon ) 
      for elec in electrons :
        t_elec = Electron(elec,'e')
        if ( t_elec.valid() ) :
          elec_list.append( t_elec )
      for jet in jets :
        t_jet = Jet(jet)
        if ( t_jet.valid() ) :
          jet_list.append(t_jet)

      ## J/psi's daughther type.
      for secVtx in secVtxs :
        if secVtx is None :
          continue
        if ( abs(secVtx.pdgId()) == 11 ) :
          #t_sV = Jpsi(secVtx, electrons)
          t_sV = Jpsi(secVtx)
          jpsi_list.append( t_sV )
        elif ( abs(secVtx.pdgId()) == 13 ) :
          #t_sV = Jpsi(secVtx, muons)
          t_sV = Jpsi(secVtx)
          jpsi_list.append( t_sV ) 
        else :
          print "error."

      # Event Selection.
      lep_list = muon_list + elec_list
      # If num. of leptons < 2, it is not dilepton events.
      if ( len(lep_list) <2 ) : continue

      # Cut Flow
      pass_flag={}
      for mode in c_modes :
        pass_flag[mode] =[False,False,False,False,False]
      for lep_idx, iso_lep1 in enumerate(lep_list) :
        for iso_lep2 in lep_list[lep_idx:] :
          mode =""
          if ( iso_lep1.type=='m' and iso_lep2.type=='m') :
            mode="MuMu"
          elif ( iso_lep1.type=='e' and iso_lep2.type=='e') :
            mode="ElEl"     
          elif ( iso_lep1.type=='m' and iso_lep2.type=='e') :
            mode="MuEl"
          
          #To make required POGs.
          ### Isolated lepton removed Jet.
          cleaned_jet_list =[]
          JpsiWithDR_5 = []
          JpsiWithvProb_dPV = []
          for jet in jet_list :
            if ( jet.jetCleaning( iso_lep1, iso_lep2) ) :
              cleaned_jet_list.append( jet )
          ### Jpsi with dR< 0.5
          for jpsi in jpsi_list :
            jpsi.JetDR(cleaned_jet_list)
            if ( jpsi.minDR< 0.5) :
              JpsiWithDR_5.append(jpsi)
            if ( jpsi.vProb>0.01 and (jpsi.l3D>0.02 and jpsi.l3D<2) ) :
              JpsiWithvProb_dPV.append( jpsi)
               
          #step1
          dilepton_mass = (iso_lep1+iso_lep2).M()
          if ( dilepton_mass<20. or iso_lep1.q* iso_lep2.q != -1 ) :
            continue
          pass_flag[mode][0] = True
          #step2
          if (dilepton_mass >76 and dilepton_mass < 106 and mode !="MuEl") :
            continue
          pass_flag[mode][1] = True
          #step3
          if ( len( cleaned_jet_list) <1 ) :
            continue
          pass_flag[mode][2] = True
          #step4
          if ( len(JpsiWithDR_5)<1 ) :
            continue
          pass_flag[mode][3] = True
          if ( len(JpsiWithvProb_dPV) <1 ) :
            continue
          pass_flag[mode][4] = True

      for mode in c_modes :
        #print pass_flag[mode]
        for idx,flag in enumerate(pass_flag[mode]) :
          if ( flag ) : 
            self.FillHist( mode, idx+1, muon_list,elec_list,jet_list,jpsi_list)

    outfile.Write()

  def Print(self) :
    print self.inputfile, self.outfilename

if __name__ == "__main__" :
  infile = sys.argv[1]
  outfile = sys.argv[2]
  ana = JpsiAna( infile, outfile )
  ana.Ana()
