from ROOT import *
from CATJpsi.JPsiAnalysis.Lepton import Lepton
class Electron(Lepton) :
  def __init__( self, Electron, type) :
    Lepton.__init__(self, Electron, type)
    self.passCV = Electron.passConversionVeto()
    self.scEta = Electron.scEta()
    self.mva = Electron.electronID('mvaTrigV0')
    self.isPF = Electron.isPF()
  def valid( self ) :
    if ( self.Pt() >20 and fabs( self.Eta()) < 2.5 and self.relIso< 0.15 and self.passCV > 0.999 and self.mva >0.5 and self.isPF>0.999 ) :
      return True
    else :
      return False 
  
