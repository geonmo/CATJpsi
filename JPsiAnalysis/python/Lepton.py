#!/usr/bin/env python
from ROOT import *
class Lepton(TLorentzVector) :
    def __init__( self, pt, eta, phi, m, q, relIso, type) :
      TLorentzVector.__init__(self)
      self.SetPtEtaPhiM( pt, eta, phi, m)
      self.q = q
      self.relIso = relIso
      self.type = type
