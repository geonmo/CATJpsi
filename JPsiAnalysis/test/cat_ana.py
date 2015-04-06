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

c_modes = ['MuMu','MuEl','ElEl']
steps = range(6)[1:]

class JpsiAna :
  inputfile = None
  hist={}
  def __init__( self, infile, outfile ) :
    self.inputfile = infile
    self.outfilename = outfile

  def BookingHist( self,output ) :
    for mode in c_modes :
      for step in steps :
        dir = "%s_s%d"%(mode,step)
        output.cd(dir)
        self.hist[(mode,step, 'muon_pt')  ] = TH1F ("muon_pt", "muon_pt", 100,0,200)
        self.hist[(mode,step, 'muon_eta') ] = TH1F ("muon_eta", "muon_eta", 100,-5,5)
        self.hist[(mode,step, 'muon_phi') ] = TH1F ("muon_phi", "muon_phi", 100,-math.pi,math.pi)
        self.hist[(mode,step, 'muon_mass')] = TH1F ("muon_mass", "muon_mass;Muon mass(GeV/c^2);Entries", 100,0.1,0.11)
        self.hist[(mode,step, 'muon_relIso')] = TH1F("muon relIso","muon_relIso",100,0,1);


        self.hist[(mode,step, 'elec_pt')  ] = TH1F ("elec_pt", "elec_pt", 100,0,200)
        self.hist[(mode,step, 'elec_eta') ] = TH1F ("elec_eta", "elec_eta", 100,-5,5)
        self.hist[(mode,step, 'elec_phi') ] = TH1F ("elec_phi", "elec_phi", 100,-math.pi,math.pi)
        self.hist[(mode,step, 'elec_mass') ] = TH1F ("elec_mass", "elec_mass;Muon mass(GeV/c^2);Entries", 100,0.1,0.11)
        self.hist[(mode,step, 'elec_relIso')] = TH1F("elec relIso","elec_relIso",100,0,1);

        self.hist[(mode,step, 'jpsi_pt')  ] = TH1F ("jpsi_pt", "jpsi_pt", 100,0,200)
        self.hist[(mode,step, 'jpsi_eta') ] = TH1F ("jpsi_eta", "jpsi_eta", 100,-5,5)
        self.hist[(mode,step, 'jpsi_phi') ] = TH1F ("jpsi_phi", "jpsi_phi", 100,-math.pi,math.pi)
        self.hist[(mode,step, 'jpsi_mass') ] = TH1F ("jpsi_mass", "jpsi_mass;Muon mass(GeV/c^2);Entries", 100,2.8,3.4)
        self.hist[(mode,step, 'jpsi_vProb') ] = TH1F ("jpsi_vProb", "jpsi_vProb;Probablity;Entries", 100,0,1)
        self.hist[(mode,step, 'jpsi_l3D') ]   = TH1F ("jpsi_l3D", "jpsi_l3D;l3D;Entries", 500,0,5)
        self.hist[(mode,step, 'l1jpsi_pt')]   = TH1F ("l1jpsi_pt", "l1jpsi_pt", 100,0,200)                                   
        self.hist[(mode,step, 'l1jpsi_eta')]  = TH1F ("l1jpsi_eta", "l1jpsi_eta", 100,-5,5)                                 
        self.hist[(mode,step, 'l1jpsi_phi')]  = TH1F ("l1jpsi_phi", "l1jpsi_phi", 100,-math.pi,math.pi)                     
        self.hist[(mode,step, 'l1jpsi_m')]    = TH1F ("l1jpsi_mass", "l1jpsi_mass;Muon mass(GeV/c^2);Entries", 100,0,100)
        self.hist[(mode,step, 'l2jpsi_pt')]   = TH1F ("l2jpsi_pt", "l2jpsi_pt", 100,0,200)                                   
        self.hist[(mode,step, 'l2jpsi_eta')]  = TH1F ("l2jpsi_eta", "l2jpsi_eta", 100,-5,5)                                 
        self.hist[(mode,step, 'l2jpsi_phi')]  = TH1F ("l2jpsi_phi", "l2jpsi_phi", 100,-math.pi,math.pi)                     
        self.hist[(mode,step, 'l2jpsi_m')]    = TH1F ("l2jpsi_mass", "l2jpsi_mass;Muon mass(GeV/c^2);Entries", 100,0,100)
        self.hist[(mode,step, 'nEvent') ] = TH1F ("nEvent", "number of Events;# of Event;Entries", 1,0,1)

  def FillHist( self, mode,step, muons,elecs,jets,jpsis ) :
      for muon in muons :
        self.hist[(mode,step,"muon_pt")].Fill( muon.Pt())
        self.hist[(mode,step,"muon_eta")].Fill( muon.Eta() )
        self.hist[(mode,step,"muon_phi")].Fill( muon.Phi() )
        self.hist[(mode,step,"muon_relIso")].Fill( muon.relIso )
      for elec in elecs :
        self.hist[(mode,step,"elec_pt")].Fill( elec.Pt())
        self.hist[(mode,step,"elec_eta")].Fill( elec.Eta() )
        self.hist[(mode,step,"elec_phi")].Fill( elec.Phi() )
        self.hist[(mode,step,"elec_mass")].Fill( elec.M() )
        self.hist[(mode,step,"elec_relIso")].Fill( elec.relIso )
      for jpsi in jpsis :
        self.hist[(mode,step, 'jpsi_pt')  ].Fill( jpsi.Pt())
        self.hist[(mode,step, 'jpsi_eta') ].Fill( jpsi.Eta())
        self.hist[(mode,step, 'jpsi_phi') ].Fill( jpsi.Phi())
        self.hist[(mode,step, 'jpsi_mass') ].Fill( jpsi.M() )
        self.hist[(mode,step, 'jpsi_vProb') ].Fill( jpsi.vProb )
        self.hist[(mode,step, 'jpsi_l3D') ].Fill(jpsi.l3D)
        l1jpsi,l2jpsi,lep1,lep2 = None, None, None, None
        if ( mode =="MuMu") :
          lep1 = muons[0]
          lep2 = muons[1]
        elif ( mode =="MuEl") :
          lep1 = muons[0]
          lep2 = elecs[0]
        elif ( mode =="ElEl") :
          lep1 = elecs[0]
          lep2 = elecs[1]
        lep1_angle = jpsi.Angle( lep1.Vect())
        lep2_angle = jpsi.Angle( lep2.Vect())
        if ( lep1_angle< lep2_angle) :
          l1jpsi = jpsi+lep1
          l2jpsi = jpsi+lep2
        else :
          l1jpsi = jpsi+lep2
          l2jpsi = jpsi+lep1
        self.hist[(mode,step, 'l1jpsi_pt')].Fill( l1jpsi.Pt())
        self.hist[(mode,step, 'l1jpsi_eta')].Fill( l1jpsi.Eta())
        self.hist[(mode,step, 'l1jpsi_phi')].Fill( l1jpsi.Phi())
        self.hist[(mode,step, 'l1jpsi_m')].Fill( l1jpsi.M())
        self.hist[(mode,step, 'l2jpsi_pt')].Fill( l2jpsi.Pt())
        self.hist[(mode,step, 'l2jpsi_eta')].Fill( l2jpsi.Eta())
        self.hist[(mode,step, 'l2jpsi_phi')].Fill( l2jpsi.Phi())
        self.hist[(mode,step, 'l2jpsi_m')].Fill( l2jpsi.M())
      self.hist[(mode,step,'nEvent')].Fill(0)
    
  def Ana( self ) :
    file = TFile(self.inputfile)
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
    catElectronLabel, catElectron = "catElectrons", Handle("std::vector<cat::Electron>")
    catSecVertexLabel, catSecVtx = "catSecVertexs", Handle("std::vector<cat::SecVertex>")
    catJetLabel, catJet = "catJets", Handle("std::vector<cat::Jet>")
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
