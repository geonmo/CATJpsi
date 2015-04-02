#!/usr/bin/env python
from ROOT import *
class Jet(TLorentzVector) :
    def __init__( self, pt, eta, phi, m, bTagCSV,dn,res_dn) :
      TLorentzVector.__init__(self)
      self.SetPtEtaPhiM( pt, eta, phi, m)
      self.bTagCSV = bTagCSV
      self.dn = dn
      self.res_dn = res_dn
    def valid( self) :
      if ( self.Pt() > 30 and fabs(self.Eta()) < 2.5 ) :
            return True
      else : 
        return False
    def isBTag( self ) :
      if ( self.bTagCSV > 0.244) :
        return True
      else :
        return False
