#!/usr/bin/env python
from ROOT import *
class Jpsi(TLorentzVector) :
  def __init__( self, jpsi) :
    TLorentzVector.__init__(self)
    self.SetPtEtaPhiM( jpsi.pt(), jpsi.eta(), jpsi.phi(), jpsi.mass())
    self.vProb_ = jpsi.vProb()
    self.l3D_ = jpsi.l3D()
    self.lxy_ = jpsi.lxy()
    self.dca_ = jpsi.dca()
    self.muonID_ = jpsi.muID()
    self.trackQuality_ = jpsi.trackQuality()
    self.sigmalxy_ = jpsi.sigmalxy()
    self.cxPtHypot_ = jpsi.cxPtHypot()
    self.cxPtAbs_ = jpsi.cxPtAbs()
    self.minDR_ = 9999.
    self.minBDR_ = 9999.

  def JetDR ( self, jet_list ) :
    for jet in jet_list :
      deltaR = self.DeltaR( jet)
      if ( deltaR < self.minDR ) : self.minDR = deltaR
  def BJetDR( self, jet_list ) :
    for jet in jet_list :
      if ( jet.isBTag()) : 
        deltaR = self.DeltaR( jet)
        if ( deltaR < self.minBDR ) : self.minBDR = deltaR
  def pt( self ) :
    return self.Pt()        
  def eta(self ) :
    return self.Eta()
  def phi(self ) :
    return self.Phi()
  def mass(self) : 
    return self.M()
  def vProb(self ) :
    return self.vProb_
  def l3D(self) :
    return self.l3D_
  def lxy(self) :
    return self.lxy_
  def dca(self) :
    return self.dca_
  def muID(self) :
    return self.muonID_
  def trackQuality(self) :
    return self.trackQuality_
  def minDR(self) :
    return  self.minDR_
  def minBDR(self) :
    return self.minBDR_
  def sigmalxy(self) :
    return self.sigmalxy_
  def cxPtHypot(self) :
    return self.cxPtHypot_
  def cxPtAbs(self ) :
    return self.cxPtAbs_

  def valid( self) :
    if (self.vProb_ > 0.01 and self.dca_<0.05 and self.l3D_>0.2) : return True
    else : return False
