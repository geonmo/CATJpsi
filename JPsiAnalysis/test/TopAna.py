#!/usr/bin/env python
from ROOT import *
gSystem.Load("libFWCoreFWLite.so");
gSystem.Load("libDataFormatsFWLite.so");
AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,math,sys
from AbsAna import AbsAna 

class TopAna(AbsAna) : 
  def __init__(self,infile, outfile) :
    AbsAna.__init__(self,infile,outfile)
  def isFromB( self, particle ) :
    if (particle is None ) :
      return False
    try :
      if ( abs(particle.pdgId()) == 5 ) :
        return True
    except :
      print "Particle is not come from 'genParticle'. Can not get pdgId()"
      return False
    flag = False
    nMother = particle.numberOfMothers()
    for idx in range(nMother) :
      if ( self.isFromB( particle.mother(idx) ) ) :
        flag = True
    return flag

  def isFromTop( self, particle ) :
    if (particle is None ) :
      return False
    try :
      if ( particle.pdgId() == 6 ) :
        return True
    except :
      print "Particle is not come from 'genParticle'. Can not get pdgId()"
      return False
    flag = False
    nMother = particle.numberOfMothers()
    for idx in range(nMother) :
      if ( self.isFromTop( particle.mother(idx) ) ) :
        flag = True
    return flag

  def isFromTopBar( self, particle ) :
    if (particle is None ) :
      return False
    try :
      if ( particle.pdgId() == -6 ) :
        return True
    except :
      print "Particle is not come from 'genParticle'. Can not get pdgId()"
      return False
    flag = False
    nMother = particle.numberOfMothers()
    for idx in range(nMother) :
      if ( self.isFromTopBar( particle.mother(idx) ) ) :
        flag = True
    return flag
