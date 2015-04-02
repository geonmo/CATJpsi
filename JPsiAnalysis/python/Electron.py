from ROOT import *
from CATJpsi.JPsiAnalysis.Lepton import Lepton
class Electron(Lepton) :
  def __init__( self, pt, eta,phi,m,q,relIso, type, mva, passConversionVeto, scEta, isPF) :
    Lepton.__init__(self, pt, eta, phi, m, q, relIso, type)
    self.passCV = passConversionVeto
    self.scEta = scEta 
    self.mva = mva
    self.isPF = isPF
  def valid( self ) :
    if ( self.Pt() >20 and fabs( self.Eta()) < 2.5 and self.relIso< 0.15 and self.passCV > 0.999 and self.mva >0.5 and self.isPF>0.999 ) :
      return True
    else :
      return False 
  
