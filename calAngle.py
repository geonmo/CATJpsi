#!/usr/bin/env python
from ROOT import *
gSystem.Load("libFWCoreFWLite.so");
gSystem.Load("libDataFormatsFWLite.so");
AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,math,sys

class Gen :
  outputfile = None
  def __init__(self,outputfile) :
    self.outputfile = outputfile
    pass

  def isFromB( self, particle ) :
    if (particle is None ) :
      return False
    if ( abs(particle.pdgId()) == 5 ) :
      return True
    flag = False
    nMother = particle.numberOfMothers()
    for idx in range(nMother) :
      if ( self.isFromB( particle.mother(idx) ) ) :
        flag = True
    return flag

  def isFromTop( self, particle ) :
    if (particle is None ) :
      return False
    if ( (particle.pdgId()) == 6 ) :
      return True
    flag = False
    nMother = particle.numberOfMothers()
    for idx in range(nMother) :
      if ( self.isFromTop( particle.mother(idx) ) ) :
        flag = True
    return flag

  def isFromTopBar( self, particle ) :
    if (particle is None ) :
      return False
    if ( (particle.pdgId()) == -6 ) :
      return True
    flag = False
    nMother = particle.numberOfMothers()
    for idx in range(nMother) :
      if ( self.isFromTopBar( particle.mother(idx) ) ) :
        flag = True
    return flag



  def isEqual( self, jpsi1, jpsi2, exact=True ) :
    rel_pt = abs((jpsi1.pt()-jpsi2.pt())/jpsi1.pt())
    d_eta = jpsi1.eta()-jpsi2.eta()
    d_phi = TVector2.Phi_mpi_pi(jpsi1.phi()-jpsi2.phi())
    dRval = sqrt( d_eta*d_eta + d_phi*d_phi)
    if ( exact ) :
      if ( rel_pt < 1e-3 and dRval < 1e-3 ) :
        return True 
    else :
      if ( rel_pt < 0.05 and dRval < 0.05 ) :
        return True
    return False
   
  def cleaning( self, jpsis,exact=True) :
    c_jpsis =[]
    for jpsi in jpsis :
      duplicate = False
      for c_jpsi in c_jpsis :
        if ( self.isEqual (jpsi,c_jpsi,exact)) :
          duplicate = True
      if ( not duplicate ) :
        c_jpsis.append(jpsi)
    return c_jpsis

  def FillNtuple(self, jpsi_tuple, jpsiList) :
    for jpsi in jpsiList :
      jpsi_tuple.Fill( jpsi.pt(), jpsi.eta(), jpsi.phi(),jpsi.mass(),0.0,0.0,jpsi.lxy(),jpsi.l3D(),jpsi.vProb(),jpsi.sigmalxy(),jpsi.dca(),jpsi.cxPtHypot(),jpsi.cxPtAbs());

  def Ana(self,infile) :
    print "inner Ana : %s"%(infile)
    file = TFile.Open(infile)
    print "Load file : "+file.GetName()
    if ( file.IsZombie()) :
      print "%s is corruct!"%(infile)
      return
    file.Close()

    events = Events(infile)
    catGenLabel, genParticle = "genParticles", Handle("std::vector<reco::GenParticle>")

    output = TFile(self.outputfile,"RECREATE")

    jpsi_info = TNtuple("jpsi_info","Information of di-muon and J/#psi pT","run:evt:dR:JPsipT:fromB:fromTop:muon1_pt:muon2_pt")
    
    for iev,event in enumerate(events):
      event.getByLabel(catGenLabel, genParticle)

      try :      
        gens = genParticle.product()
      except RuntimeError :
        print "Warning! No genParticle is found."
        continue

      gen_jpsis =[]
      gen_muons =[]
      ori_gen_jpsis=[]
      c_gen_jpsis=[]
      evt =0
      run =0

      for gen in gens :
        if ( gen.pdgId() == 443 and gen.pt()>1.0 and abs(gen.eta())<2.5  and gen.mass()>3 and gen.mass()<3.2) :
          if ( gen.numberOfDaughters() ==2 ) :
            if ( abs(gen.daughter(0).pdgId())==13 and abs(gen.daughter(1).pdgId())==13 ) :
              gen_jpsis.append( gen )
              #print "GEN Found"
        if ( abs(gen.pdgId()) == 13 and  gen.pt()>1 and abs(gen.eta())<2.5 and gen.status() ==1 ) :
          gen_muons.append( gen)

      c_gen_jpsis = self.cleaning( gen_jpsis,True)
      if ( len(c_gen_jpsis) < 1 or len(gen_muons)< 2 ) :
        continue


      if ( len(c_gen_jpsis)>0 and len(gen_muons)> 1 ) :
        evt =  event.eventAuxiliary().id().event() 
        run =  event.eventAuxiliary().id().run()
        print "Run : %d , Event : %d"%(run,evt) 
      else :
        continue

      for genJpsi in c_gen_jpsis :
        first_daughter  = TLorentzVector()
        first_daughter.SetPtEtaPhiM( genJpsi.daughter(0).pt(), genJpsi.daughter(0).eta(), genJpsi.daughter(0).phi(), genJpsi.daughter(0).mass())
        second_daughter = TLorentzVector()
        second_daughter.SetPtEtaPhiM( genJpsi.daughter(1).pt(), genJpsi.daughter(1).eta(), genJpsi.daughter(1).phi(), genJpsi.daughter(1).mass())
        dR = ( first_daughter.DeltaR( second_daughter))
        jpsi_isFromB = self.isFromB(genJpsi)
        jpsi_isFromTop = self.isFromTop(genJpsi)
        jpsi_isFromTopBar = self.isFromTopBar(genJpsi)
     
        iso = 0 
        muon1_pt =0
        muon2_pt =0
        for genMuon in gen_muons :
          muon_isFromTop = self.isFromTop(genMuon)
          muon_isFromTopBar = self.isFromTopBar(genMuon)
          if ( muon_isFromTop or muon_isFromTopBar ) :
            iso = iso+1
            if ( (jpsi_isFromTop and muon_isFromTop) or (jpsi_isFromTopBar and muon_isFromTopBar) ) :
              muon1_pt = genMuon.pt()
            else :
              muon2_pt = genMuon.pt()
            if ( iso == 2 ) : 
              jpsi_info.Fill( run, evt, dR, genJpsi.pt(), jpsi_isFromB, (jpsi_isFromTop or jpsi_isFromTopBar), muon1_pt, muon2_pt)
              break
        

    output.Write()
    output.Close()
    
      
   
if __name__ == "__main__" :
  gen = Gen(sys.argv[2])
  gen.Ana(sys.argv[1])
