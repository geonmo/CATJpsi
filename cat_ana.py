#!/usr/bin/env python
from ROOT import *
from Lepton import Lepton
from Muon import Muon
from Jet import Jet
from Electron import Electron

gSystem.Load("libFWCoreFWLite.so");
gSystem.Load("libDataFormatsFWLite.so");
AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,math
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
        #self.muon_mass = TH1F ("muon_mass", "muon_mass;Muon mass(GeV/c^2);Entries", 100,0.1,0.11)
        self.hist[(mode,step, 'muon_relIso')] = TH1F("muon relIso","muon_relIso",100,0,1);
        """
        self.dimuon_pt[mode] = TH1F ("dimuon_pt", "dimuon pt", 100,0,200)
        self.dimuon_eta[mode] = TH1F ("dimuon_eta", "dimuon eta", 100,-5,5)
        self.dimuon_phi[mode] = TH1F ("dimuon_phi", "dimuon phi", 100,-math.pi,math.pi)
        self.dimuon_mass[mode] = TH1F ("dimuon_mass", "dimuon mass;Muon mass(GeV/c^2);Entries", 100,0,5)
        self.num_jet[mode] = TH1F ("Number of Jets", "num_jet", 10,0,10)
        """

  def FillHist( self, modes,step, muons,elecs,jets ) :
    for mode in modes :
        for muon in muons :
          self.hist[(mode,step,"muon_pt")].Fill( muon.Pt())
          self.hist[(mode,step,"muon_eta")].Fill( muon.Eta() )
          self.hist[(mode,step,"muon_phi")].Fill( muon.Phi() )
          self.hist[(mode,step,"muon_relIso")].Fill( muon.relIso )
    
  def Ana( self ) :
    file = TFile(self.inputfile)
    print "Load file : "+file.GetName()
    if ( file.IsZombie()) :
      print "%s is corruct!"%(self.inputfile)
      return
    file.Close()

    outfile = TFile("hist/"+self.outfilename+".root","RECREATE")
    for mode in c_modes :
      #outfile.mkdir(mode)
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
        print "Warning! Among of one's product is null."
      if (muons.size()>2 or electrons.size()>2) :
        secVtxs = catSecVtx.product()
      else :
        secVtxs = []

      modes =[]

      muon_list =[]
      elec_list =[]
      jet_list =[]
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


      if len( muon_list ) >= 2 :
        modes.append("MuMu")
      if (len(muon_list) >=1 and len(elec_list ) >= 1) :
        modes.append("MuEl")
      if len(elec_list) >=2 : 
        modes.append("ElEl")

      # Cut Flow
      pass_step = 0      
      #self.FillHist( modes, 1, isolated_muon )
      lep_list = muon_list + elec_list
      if ( len(lep_list) <2 ) : continue
      f_lep = lep_list[0]
      for lep in lep_list[1:] :
        #step1
        if ( (f_lep + lep).M() > 20. and f_lep.q* lep.q == -1 ) :
          pass_step = 1
          #step2 
          if  (f_lep+lep).M() <76. or (f_lep+lep).M() > 106.  :
            pass_step = 2
            #step3
            cleaned_jet_list =[]
            for jet in jet_list :
              if ( jet.jetCleaning( f_lep, lep) ) :
                cleaned_jet_list.append( jet )
                #print "jet cleaned"
            if ( len( cleaned_jet_list) >=1 ) :
              pass_step = 3
              #step4
              if ( len( secVtxs ) > 0 ) :
                pass_step = 4

      for ps in range(pass_step+1)[1:] : 
        self.FillHist( modes, ps, muon_list,elec_list,jet_list)

    outfile.Write()

  def Print(self) :
    print self.inputfile, self.outfilename

if __name__ == "__main__" :
  #signal_path = "/cms/data/xrd/store/user/geonmo/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola/cat-ttbar_Nov20_2014/141120_114344/0000/"
  signal_path = "./"
  #files = os.listdir(signal_path)
  files = ["catTuple_1.root"]
  for i,file in enumerate(files) :
    ana = JpsiAna( signal_path+file, "hist_%04d"%(i) )
    ana.Ana()
    ana.Print()
    #ana.__del__()
