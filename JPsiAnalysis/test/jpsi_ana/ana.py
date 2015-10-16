#!/usr/bin/env python
from template.pyROOT_TreeAna import TreeAna
from ROOT import *
import os
class JpsiAna(TreeAna) :
  def __init__(self) :
    TreeAna.__init__(self)
    if  not os.access(("hist"), os.R_OK) :
      os.makedirs("hist")
    self.hist['jpsi_mass'] = TH1F("jpsi_mass","jpsi mass",100,3,3.2)
    self.hist['ljpsi1_mass'] = TH1F("ljpsi1_mass","ljpsi1 mass",85 ,10 , 200)
    self.hist['ljpsi2_mass'] = TH1F("ljpsi2_mass","ljpsi2 mass",85 ,10 , 200)
  def Ana( self, mychain) :
    if ( mychain.step !=5 ) :
      return 
    for idx in range( mychain.jpsi_mass.size() ) :
      if ( mychain.jpsi_mva[idx] > -0.0360 and mychain.jpsi_minBDR[idx] < 0.5 ) :
        self.hist['jpsi_mass'].Fill( mychain.jpsi_mass[idx], mychain.pileupWeight )
        self.hist['ljpsi1_mass'].Fill( mychain.ljpsi1_m[idx], mychain.pileupWeight)
        self.hist['ljpsi2_mass'].Fill( mychain.ljpsi2_m[idx], mychain.pileupWeight)
   
## If directly use ths script,
if __name__ == "__main__" :
  ana = JpsiAna( )
  ana.Loop()
  ana.Write()
