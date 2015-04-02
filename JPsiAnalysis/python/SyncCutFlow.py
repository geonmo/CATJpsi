#!/usr/bin/env python
from ROOT import *
class Flag :
  def __init__(self) :
    self.mm  = [False,False,False,False,False]
    self.ee  = [False,False,False,False,False]
    self.me  = [False,False,False,False,False]

class CutFlow :
  def __init__(self, lep1, lep2, flag ) :
    self.lep1 = lep1
    self.lep2 = lep2
    self.zBoson = lep1+lep2
    self.flag = flag
    self.nbJets = 0

  def SetFlag( self, idx) :
    if   ( self.lep1.type == 'm' and self.lep2.type == 'm' ) :
      self.flag.mm[idx] = True
    elif ( self.lep1.type == 'm' and self.lep2.type == 'e' ) :
      self.flag.me[idx] = True
    elif ( self.lep1.type == 'e' and self.lep2.type == 'm' ) :
      self.flag.me[idx] = True
    elif ( self.lep1.type == 'e' and self.lep2.type == 'e' ) :
      self.flag.ee[idx] = True
  def JetCleaning( self, jet_list, lep_list) :
    cleaned_jets =[]
    for jet in jet_list :
      jet_keep_flag = True
      for lep in lep_list : 
        if ( jet.DeltaR( lep ) < 0.5 ) :
          jet_keep_flag = False
      if( jet_keep_flag ) :
          cleaned_jets.append( jet ) 
    return cleaned_jets
    
  def Print(self) :
    print self.lep1.M(),"  ",self.lep2.M(),"  ",self.zBoson.M()
  def Step1( self ) :
    if ( self.lep1.valid() and self.lep2.valid() and self.zBoson.M() >20 and int(self.lep1.q*self.lep2.q) == -1 ) :
      self.SetFlag(0)
      return True
    else :
      return False
  def Step2( self ) :
    if ( self.lep1.type == self.lep2.type ) :
      if ( self.zBoson.M() < 76 or self.zBoson.M() > 106 ) :
        self.SetFlag(1)
        return True
      else :
        return False
    else :
      self.SetFlag(1)
      return True
  def Step3( self, nJets ) :
    if ( nJets >=2 ) :
      self.SetFlag(2)
      return True
    else :
      return False
  def Step4( self, met_pt) :
    if ( self.lep1.type == self.lep2.type ) :
      if( met_pt > 40 ) :
        self.SetFlag(3)
        return True
      else :
        return False
    else :
      self.SetFlag(3)
      return True
  def Step5( self, jets_list ) :
    for jet in jets_list :
      if ( jet.isBTag() ) :
        self.nbJets = self.nbJets + 1
    if ( self.nbJets >=1 ) : 
      self.SetFlag(4)
      return True
    else :
      return False
