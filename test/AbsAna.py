#!/usr/bin/env python
from ROOT import *
gSystem.Load("libFWCoreFWLite.so");
gSystem.Load("libDataFormatsFWLite.so");
AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,math,sys,gc

class AbsAna :
  outfile = None
  infile = None
  def __init__(self,infile, outfile) :
    self.infile = infile
    self.outfile = outfile
    pass

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

  def isEqual( self, particle1, particle2, exact=True ) :
    rel_pt = abs((particle1.pt()-particle2.pt())/particle1.pt())
    d_eta = particle1.eta()-particle2.eta()
    d_phi = TVector2.Phi_mpi_pi(particle1.phi()-particle2.phi())
    dRval = sqrt( d_eta*d_eta + d_phi*d_phi)
    if ( exact ) :
      if ( rel_pt < 1e-3 and dRval < 1e-3 ) :
        return True 
    else :
      if ( rel_pt < 0.05 and dRval < 0.05 ) :
        return True
    return False
   
  def cleaning( self, particles,exact=True) :
    c_particles =[]
    for particle in particles :
      duplicate = False
      for c_particle in c_particles :
        if ( self.isEqual (particle,c_particle,exact)) :
          duplicate = True
      if ( not duplicate ) :
        c_particles.append(particle)
    return c_particles

  def checkInfile(self) :
    print "Checking input file : %s"%(self.infile)
    try :
      file = TFile.Open(self.infile)
    except :
      print "File open error!"
      return False

    if ( file.IsZombie()) :
      print "%s is corruct!"%(infile)
      return False
    print "Open is success. It is a right root file."
    file.Close()
    return True
    

  def Ana(self) :
    if ( not self.checkInfile() ) :
      sys.exit(-1)
    output = TFile(self.outfile, "RECREATE")
    output.Close()
    return
   
if __name__ == "__main__" :
  infile = sys.argv[1]
  outfile = sys.argv[2]
  gen = AbsAna(infile, outfile)
  gen.Ana()
  gc.collect()
