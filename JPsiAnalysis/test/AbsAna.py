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
  events = None
  def __init__(self,infiles, outfile) :
    self.infiles = self.parsingInfile(infiles)
    self.outfile = outfile
    if ( not self.checkInfile() ) :
      print "ROOT file is not normal. Terminate analysis code."
      sys.exit(-1)
    pass
    self.events = Events(self.infiles)

  def parsingInfile(self, infiles) :
    infiles = infiles.strip()
    infile_list = infiles.replace("\"","").replace("\'","").replace("\n","").split(",")
    return infile_list
  def isEqual( self, particle1, particle2, exact=True ) :
    rel_pt = abs((particle1.pt()-particle2.pt())/particle1.pt())
    d_eta = particle1.eta()-particle2.eta()
    d_phi = TVector2.Phi_mpi_pi(particle1.phi()-particle2.phi())
    dRval = sqrt( d_eta*d_eta + d_phi*d_phi)
    if ( exact ) :
      if ( rel_pt < 1e-3 and dRval < 1e-3 ) :
        return True 
    else :
      if ( rel_pt < 0.05 and dRval < 0.15 ) :
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
    print self.infiles
    for infile in self.infiles :
      print "Checking input file : %s"%(infile)
      try :
        file = TFile.Open(infile)
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
    pass
   
if __name__ == "__main__" :
  infile = sys.argv[1]
  outfile = sys.argv[2]
  gen = AbsAna(infile, outfile)
  gen.Ana()
  gc.collect()
