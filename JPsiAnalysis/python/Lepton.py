#!/usr/bin/env python
from ROOT import *
class Lepton(TLorentzVector) :
    def __init__( self, lepton, type ) :
      TLorentzVector.__init__(self)
      self.SetPtEtaPhiM( lepton.pt(), lepton.eta(), lepton.phi(), lepton.mass() )
      self.q = lepton.charge()
      self.relIso = lepton.relIso()
      self.type = type
