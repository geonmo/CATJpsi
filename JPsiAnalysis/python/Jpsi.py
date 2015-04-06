#!/usr/bin/env python
from ROOT import *
class Jpsi(TLorentzVector) :
  def __init__( self, jpsi) :
    TLorentzVector.__init__(self)
    self.SetPtEtaPhiM( jpsi.pt(), jpsi.eta(), jpsi.phi(), jpsi.mass())
    self.vProb = jpsi.vProb()
    self.l3D = jpsi.l3D()
    self.lxy = jpsi.lxy()
    self.dca = jpsi.dca()
    self.muonID = jpsi.muID()
    self.trackIndex = jpsi.trackIndex()
    self.cxPtHypot = jpsi.cxPtHypot()
    self.minDR = 9999.

  def JetDR ( self, jet_list ) :
    for jet in jet_list :
      deltaR = self.DeltaR( jet)
      if ( deltaR < self.minDR ) : self.minDR = deltaR
