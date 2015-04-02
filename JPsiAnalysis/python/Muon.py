from ROOT import *
from CATJpsi.JPsiAnalysis.Lepton import Lepton
class Muon(Lepton) :
  def __init__( self, pt, eta,phi,m,q,relIso, type, isLoose) :
    Lepton.__init__(self, pt, eta, phi, m, q, relIso, type)
    self.isLoose = isLoose
  def valid( self ) :
    if ( self.Pt() >20 and fabs( self.Eta()) < 2.4 and self.relIso< 0.15 and self.isLoose > 0.999 ) :
      return True
    else :
      return False 
