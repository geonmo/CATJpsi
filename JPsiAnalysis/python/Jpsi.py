#!/usr/bin/env python
from ROOT import *
class Jpsi(TLorentzVector) :
  def __init__( self, pt, eta, phi, m, vProb, l3D ) :
    TLorentzVector.__init__(self)
    self.SetPtEtaPhiM( pt, eta, phi, m)
    self.vProb = vProb
    self.l3D = l3D
    self.minDR = 9999.

  def JetDR ( self, jet_list ) :
    for jet in jet_list :
      deltaR = self.DeltaR( jet)
      if ( deltaR < self.minDR ) : self.minDR = deltaR
    #return deltaR
