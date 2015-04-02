#!/usr/bin/env python
from ROOT import *
from CATJpsi.JPsiAnalysis.Muon import Muon
from CATJpsi.JPsiAnalysis.Electron import Electron
from CATJpsi.JPsiAnalysis.Lepton import Lepton
from CATJpsi.JPsiAnalysis.Jet import Jet
from CATJpsi.JPsiAnalysis.CutFlow import *
from CATJpsi.JPsiAnalysis.Jpsi import Jpsi

class EventContent :
  def __init__(self, ev,mode) :
    self.muon_list =[]
    self.elec_list =[]
    self.jet_list = []
    self.jpsimm_list =[]
    self.jpsiee_list =[]
    self.cleanedJet_list=[]
    self.lep_list=[]
    self.jpsi_list =[]
    self.ev = ev 
    self.mode = mode
    self.InitLepton()
    self.InitJet()
    self.InitJpsi()
    if len(self.lep_list) >=2 :
      self.lep1 = self.lep_list[0]
      self.lep2 = self.lep_list[1]
    else :
      pass

  def InitMuon( self ) :
    for i in range( len( self.ev.muons_pt) ) :
      temp_muon = Muon( self.ev.muons_pt[i], self.ev.muons_eta[i], self.ev.muons_phi[i], self.ev.muons_m[i], self.ev.muons_q[i], self.ev.muons_relIso[i], 'm',self.ev.muons_isLoose[i])
      self.muon_list.append( temp_muon )
  def InitElectron( self ) :
      for i in range( len( self.ev.electrons_pt) ) :
        temp_elec =  Electron( self.ev.electrons_pt[i], self.ev.electrons_eta[i], self.ev.electrons_phi[i], self.ev.electrons_m[i], self.ev.electrons_q[i], self.ev.electrons_relIso[i], 'e', self.ev.electrons_mva[i], self.ev.electrons_conversionVeto[i], self.ev.electrons_scEta[i], self.ev.electrons_isPF[i] )
        self.elec_list.append( temp_elec)

  def InitLepton( self ) :
    if ( self.mode == "MuMu") :
      self.InitMuon()
      self.lep_list = self.muon_list
    elif ( self.mode == "ElEl" ) :
      self.InitElectron() 
      self.lep_list = self.elec_list
    elif ( self.mode == "MuEl" ) :
      self.InitMuon()
      self.InitElectron()
      self.lep_list = self.muon_list + self.elec_list

  def InitJet( self ) :
    for i in range( len ( self.ev.jets_pt) ) :
      self.jet_list.append( Jet(self.ev.jets_pt[i], self.ev.jets_eta[i], self.ev.jets_phi[i], self.ev.jets_m[i], self.ev.jets_bTagCSV[i], self.ev.jets_dn[i], self.ev.jets_resDn[i]))
    self.JetCleaning()

  def InitJpsi( self): 
    for i in range( len( self.ev.jpsiMuMu_pt) ) :
      temp_jpsi = Jpsi( self.ev.jpsiMuMu_pt[i], self.ev.jpsiMuMu_eta[i], self.ev.jpsiMuMu_phi[i], self.ev.jpsiMuMu_m[i], self.ev.jpsiMuMu_vProb[i], self.ev.jpsiMuMu_l3D[i])
      temp_jpsi.JetDR( self.cleanedJet_list )
      self.jpsimm_list.append(temp_jpsi)

    for i in range( len( self.ev.jpsiElEl_pt) ) :
      temp_jpsi =  Jpsi( self.ev.jpsiElEl_pt[i], self.ev.jpsiElEl_eta[i], self.ev.jpsiElEl_phi[i], self.ev.jpsiElEl_m[i], self.ev.jpsiElEl_vProb[i], self.ev.jpsiElEl_l3D[i])
      temp_jpsi.JetDR( self.cleanedJet_list )
      self.jpsiee_list.append( temp_jpsi)

    self.jpsi_list = self.jpsimm_list+ self.jpsiee_list
    
  def JetCleaning( self ) :
    cleaned_jets =[]
    for jet in self.jet_list :
      jet_keep_flag = True
      for lep in self.lep_list : 
        if ( jet.DeltaR( lep ) < 0.5 ) :
          jet_keep_flag = False
      if( jet_keep_flag ) :
          cleaned_jets.append( jet )
    self.cleanedJet_list = cleaned_jets 
    
