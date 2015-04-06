#!/usr/bin/env python
from ROOT import *
class Jet(TLorentzVector) :
    def __init__( self, jet) :
      TLorentzVector.__init__(self)
      self.SetPtEtaPhiM( jet.pt(), jet.eta(), jet.phi(), jet.mass())
      self.bTagCSV = jet.bDiscriminator('combinedSecondaryVertexBJetTags')
      self.shiftEnDown = jet.shiftedEnDown()
      self.shiftEnUp   = jet.shiftedEnUp()
      self.smearedRes = jet.smearedRes()
      self.smearedResDown = jet.smearedResDown()
      self.smearedResUp = jet.smearedResUp()
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
    def jetCleaning( self, iso_lep1, iso_lep2 ) :
      if ( self.DeltaR( iso_lep1) < 0.5 or self.DeltaR(iso_lep2) < 0.5 ) :
        return True
      else : 
        return False 
