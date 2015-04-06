from ROOT import *
from CATJpsi.JPsiAnalysis.Lepton import Lepton
class Muon(Lepton) :
  def __init__( self, muon,type) :
    Lepton.__init__(self, muon,type)
    self.isLoose = muon.isLooseMuon()
    self.isTight = muon.isTightMuon()
  def valid( self ) :
    if ( self.Pt() >20 and fabs( self.Eta()) < 2.4 and self.relIso< 0.15 and self.isLoose > 0.999 ) :
      return True
    else :
      return False 
