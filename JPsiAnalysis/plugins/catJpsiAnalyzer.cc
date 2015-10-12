#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATJpsi/JPsiAnalysis/plugins/catJpsiAnalyzer.h"

#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

bool compareGenDR(const std::pair< float, cat::SecVertex* >& left, const std::pair< float, cat::SecVertex* >& right ){
  return left.first < right.first;
}


catJpsiAnalyzer::catJpsiAnalyzer(const edm::ParameterSet& iConfig)
{
  muonToken_ = iConfig.getParameter<edm::InputTag>("muons");
  elecToken_ = iConfig.getParameter<edm::InputTag>("electrons");
  jetToken_  = iConfig.getParameter<edm::InputTag>("jets");
  metToken_  = iConfig.getParameter<edm::InputTag>("mets");     
  secToken_  = iConfig.getParameter<edm::InputTag>("secVtxs");     
  vtxToken_  = iConfig.getParameter<edm::InputTag>("vertices");
  mcLabel_   = iConfig.getParameter<edm::InputTag>("mcLabel");

  sysEnergy_      = iConfig.getParameter<int >("energy");

  edm::Service<TFileService> fs;

  matchedtree_ = fs->make<TTree>("matched", "matched");
  matchedtree_->Branch("jpsi_pt",&b_mjpsi_pt, "jpsi_pt/F");
  matchedtree_->Branch("jpsi_eta",&b_mjpsi_eta,"jpsi_eta/F");
  matchedtree_->Branch("jpsi_phi",&b_mjpsi_phi,"jpsi_phi/F");
  matchedtree_->Branch("jpsi_m",&b_mjpsi_m,"jpsi_m/F");
  matchedtree_->Branch("jpsi_vProb",&b_mjpsi_vProb,"jpsi_vProb/F");
  matchedtree_->Branch("jpsi_lxy",&b_mjpsi_lxy,"jpsi_lxy/F");
  matchedtree_->Branch("jpsi_dca",&b_mjpsi_dca,"jpsi_dca/F");
  matchedtree_->Branch("jpsi_muID",&b_mjpsi_muID,"jpsi_muID/I");
  matchedtree_->Branch("jpsi_trackQuality",&b_mjpsi_trackQuality,"jpsi_trackQuality/I");
  matchedtree_->Branch("jpsi_minJetDR",&b_mjpsi_minJetDR,"jpsi_minJetDR/F");
  matchedtree_->Branch("jpsi_minBJetDR",&b_mjpsi_minBJetDR,"jpsi_minBJetDR/F");

  unmatchedtree_ = fs->make<TTree>("unmatched", "unmatched");
  unmatchedtree_->Branch("jpsi_pt",&b_ujpsi_pt,"jpsi_pt/F");
  unmatchedtree_->Branch("jpsi_eta",&b_ujpsi_eta,"jpsi_eta/F");
  unmatchedtree_->Branch("jpsi_phi",&b_ujpsi_phi,"jpsi_phi/F");
  unmatchedtree_->Branch("jpsi_m",&b_ujpsi_m,"jpsi_m/F");
  unmatchedtree_->Branch("jpsi_vProb",&b_ujpsi_vProb,"jpsi_vProb/F");
  unmatchedtree_->Branch("jpsi_lxy",&b_ujpsi_lxy,"jpsi_lxy/F");
  unmatchedtree_->Branch("jpsi_dca",&b_ujpsi_dca,"jpsi_dca/F");
  unmatchedtree_->Branch("jpsi_muID",&b_ujpsi_muID,"jpsi_muID/I");
  unmatchedtree_->Branch("jpsi_trackQuality",&b_ujpsi_trackQuality,"jpsi_trackQuality/I");
  unmatchedtree_->Branch("jpsi_minJetDR",&b_ujpsi_minJetDR,"jpsi_minJetDR/F");
  unmatchedtree_->Branch("jpsi_minBJetDR",&b_ujpsi_minBJetDR,"jpsi_minBJetDR/F");
}
catJpsiAnalyzer::~catJpsiAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//
bool catJpsiAnalyzer::isEqual( TLorentzVector& gen, TLorentzVector& reco, bool exact = false ) {
  float rel_pt = TMath::Abs( gen.Pt()-reco.Pt())/gen.Pt();
  float dR = gen.DeltaR( reco );
  if ( exact  && rel_pt < 1e-3 && dR < 1e-3 ) return true;
  else if ( !exact && rel_pt < 0.05 && dR < 0.15 ) return true;
  else return false;
}

// ------------ method called for each event  ------------
void catJpsiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  runOnMC_ = !iEvent.isRealData();

  edm::Handle<reco::VertexCollection> vertices;
  if ( ! iEvent.getByLabel(vtxToken_, vertices)) return;
  //const reco::Vertex &PV = vertices->front();
  if (vertices->empty()) return;



  edm::Handle<edm::View<cat::Muon> > muons;
  iEvent.getByLabel(muonToken_, muons);

  edm::Handle<edm::View<cat::Electron> > electrons;
  iEvent.getByLabel(elecToken_, electrons);

  edm::Handle<edm::View<cat::Jet> > jets;
  iEvent.getByLabel(jetToken_, jets);

  edm::Handle<edm::View<cat::MET> > mets;
  iEvent.getByLabel(metToken_, mets);
  edm::Handle<reco::GenParticleCollection> genParticles;

  edm::Handle<edm::View<cat::SecVertex> > svertxs ;
  iEvent.getByLabel(secToken_, svertxs);

  std::vector<TLorentzVector> genJpsis;
  if (runOnMC_){
    iEvent.getByLabel(mcLabel_,genParticles); 
    for (const reco::GenParticle & g : *genParticles){
      if (g.pdgId() == 443 && g.pt() >1.0 && TMath::Abs(g.eta())<2.5) {
        if ( g.numberOfDaughters() ==2) {
          if( g.daughter(0)->pt()>1.0 && g.daughter(1)->pt()>1.0 && TMath::Abs(g.daughter(0)->eta())<2.5 && abs(g.daughter(1)->eta())<2.5 && g.daughter(0)->status()==1 && g.daughter(1)->status()==1 ) {
            //if( isFromB( g) ) {
            genJpsis.push_back(TLorentzVector(g.px(), g.py(), g.pz(), g.energy()));
            // }
          }
        }
      }
    }
  }
  vector<cat::Muon> selectedMuons = selectMuons( muons.product() );
  vector<cat::Electron> selectedElectrons = selectElecs( electrons.product() );

  vector<TLorentzVector> recolep; 
  for (auto lep : selectedMuons){ recolep.push_back(lep.tlv()); }
  for (auto lep : selectedElectrons){ recolep.push_back(lep.tlv()); }
  if (recolep.size() != 2){
    //ttree_->Fill();
    return;
  }
  vector<cat::Jet> selectedJets = selectJets( jets.product(), recolep );
  vector<cat::Jet> selectedBJets = selectBJets( selectedJets );

  //  printf("selectedMuons %lu, selectedElectrons %lu, selectedJets %lu, selectedBJets %lu\n",selectedMuons.size(), selectedElectrons.size(), selectedJets.size(), selectedBJets.size() );

  auto svtxs = svertxs.product();
  std::vector< cat::SecVertex > new_catJpsi;
  std::vector< std::pair<float, cat::SecVertex*> > matched_catJpsi_pair;
  std::vector< cat::SecVertex > unmatched_catJpsi;

  for ( auto catJpsi = svtxs->begin(), end= svtxs->end() ; catJpsi != end ; ++catJpsi) {
    bool duplicate = false;
    for ( auto comp_catJpsi = new_catJpsi.begin(), comp_end = new_catJpsi.end() ; comp_catJpsi != comp_end ;++comp_catJpsi) {
      TLorentzVector jpsi_tlv = TLorentzVector( catJpsi->px(), catJpsi->py(), catJpsi->pz(), catJpsi->energy());
      TLorentzVector comp_jpsi_tlv = TLorentzVector( comp_catJpsi->px(), comp_catJpsi->py(), comp_catJpsi->pz(), comp_catJpsi->energy());
      if ( isEqual( jpsi_tlv, comp_jpsi_tlv) ) {
        duplicate = true ; 
      }
    }
    if ( !duplicate ) new_catJpsi.push_back( *catJpsi) ;
  }

  if ( new_catJpsi.size() >0 ) { 
    for ( auto catJpsi = new_catJpsi.begin(), end= new_catJpsi.end() ; catJpsi != end ; ++catJpsi) {
      bool matched = false;
      float min_deltaR=999.;
      TLorentzVector jpsi_tlv = TLorentzVector( catJpsi->px(), catJpsi->py(), catJpsi->pz(), catJpsi->energy());
      for ( auto& genJpsi : genJpsis) {
        if ( isEqual( jpsi_tlv, genJpsi) ) {
          matched = true; 
          float deltaR = jpsi_tlv.DeltaR( genJpsi);
          if ( min_deltaR > deltaR) min_deltaR = deltaR;
        }
      }
      if ( matched) matched_catJpsi_pair.push_back( std::make_pair<float&, cat::SecVertex*>(min_deltaR, &*catJpsi)  );
      else unmatched_catJpsi.push_back( *catJpsi);
    }
    std::sort( matched_catJpsi_pair.begin(), matched_catJpsi_pair.end(), compareGenDR);

    if ( matched_catJpsi_pair.size()>0) {
      cat::SecVertex* catJpsi = matched_catJpsi_pair[0].second;
      int mu_id = catJpsi->muID();
      float min_DR = 999.;
      float min_BDR = 999.;
      TLorentzVector jpsi_tlv = TLorentzVector( catJpsi->px(), catJpsi->py(), catJpsi->pz(), catJpsi->energy());
      for (auto jet = selectedJets.begin(), end = selectedJets.end(); jet != end; ++jet){
        TLorentzVector jet_tlv = jet->tlv();
        float deltaR = jpsi_tlv.DeltaR( jet_tlv) ;
        if ( min_DR > deltaR ) min_DR = deltaR; 
      }
      for (auto bjet = selectedBJets.begin(), end = selectedBJets.end(); bjet != end; ++bjet){
        TLorentzVector bjet_tlv = bjet->tlv();
        float deltaR = jpsi_tlv.DeltaR( bjet_tlv);
        if ( min_BDR > deltaR ) min_BDR = deltaR; 
      }
      if ( min_BDR>0 && min_BDR<10) {
        b_mjpsi_pt=(catJpsi->pt());
        b_mjpsi_eta=(catJpsi->eta());
        b_mjpsi_phi=(catJpsi->phi());
        b_mjpsi_m=(catJpsi->mass());
        b_mjpsi_vProb=(catJpsi->vProb());
        b_mjpsi_lxy=(catJpsi->lxy());
        b_mjpsi_dca=(catJpsi->dca());
        b_mjpsi_muID=( (mu_id&1) + (mu_id&2));
        b_mjpsi_trackQuality=(catJpsi->trackQuality());
        b_mjpsi_minJetDR=( min_DR ); 
        b_mjpsi_minBJetDR=( min_BDR );
        matchedtree_->Fill(); 
      }
    }

    for ( auto catJpsi = unmatched_catJpsi.begin(), end= unmatched_catJpsi.end() ; catJpsi != end ; ++catJpsi){
      int mu_id = catJpsi->muID();
      float min_DR = 999.;
      float min_BDR = 999.;
      TLorentzVector jpsi_tlv = TLorentzVector( catJpsi->px(), catJpsi->py(), catJpsi->pz(), catJpsi->energy());
      for (auto jet = selectedJets.begin(), end = selectedJets.end(); jet != end; ++jet){
        TLorentzVector jet_tlv = jet->tlv();
        float deltaR = jpsi_tlv.DeltaR( jet_tlv) ;
        if ( min_DR > deltaR ) min_DR = deltaR; 
      }
      for (auto bjet = selectedBJets.begin(), end = selectedBJets.end(); bjet != end; ++bjet){
        TLorentzVector bjet_tlv = bjet->tlv();
        float deltaR = jpsi_tlv.DeltaR( bjet_tlv);
        if ( min_BDR > deltaR ) min_BDR = deltaR; 
      }
      if ( min_BDR>0 && min_BDR<10) {
        b_ujpsi_pt=(catJpsi->pt());
        b_ujpsi_eta=(catJpsi->eta());
        b_ujpsi_phi=(catJpsi->phi());
        b_ujpsi_m=(catJpsi->mass());
        b_ujpsi_vProb=(catJpsi->vProb());
        b_ujpsi_lxy=(catJpsi->lxy());
        b_ujpsi_dca=(catJpsi->dca());
        b_ujpsi_muID=( (mu_id&1) + (mu_id&2));
        b_ujpsi_trackQuality=(catJpsi->trackQuality());
        b_ujpsi_minJetDR=( min_DR ); 
        b_ujpsi_minBJetDR=( min_BDR );
        unmatchedtree_->Fill(); 
      }
    } 
  }
}

const reco::Candidate* catJpsiAnalyzer::getLast(const reco::Candidate* p)
{
  for ( size_t i=0, n=p->numberOfDaughters(); i<n; ++i )
  {
    const reco::Candidate* dau = p->daughter(i);
    if ( p->pdgId() == dau->pdgId() ) return getLast(dau);
  }
  return p;
}

vector<cat::Muon> catJpsiAnalyzer::selectMuons(const edm::View<cat::Muon>* muons )
{
  vector<cat::Muon> selmuons;
  for (auto mu : *muons) {
    if (!mu.isTightMuon()) continue;
    if (mu.pt() <= 20.) continue;
    if (fabs(mu.eta()) >= 2.4) continue;
    if (mu.relIso(0.4) >= 0.12) continue;
    selmuons.push_back(mu);
  }

  return selmuons;
}

vector<cat::Electron> catJpsiAnalyzer::selectElecs(const edm::View<cat::Electron>* elecs )
{
  vector<cat::Electron> selelecs;
  /*
     for (auto el : *elecs) {
     if (!el.electronID("eidLoose")) continue;
     if (!el.passConversionVeto()) continue;
     if (!el.isPF()) continue;
     if (el.pt() <= 20.) continue;
     if ((fabs(el.scEta()) <= 1.4442) && (el.relIso(0.3) >= 0.1649)) continue;
     if ((fabs(el.scEta()) >= 1.566) && (el.relIso(0.3) >= 0.2075)) continue;
     if ((fabs(el.scEta()) > 1.4442) && (fabs(el.scEta()) < 1.566)) continue;
     if (fabs(el.eta()) >= 2.5) continue;
     if (el.pt() < 5) continue;
  //printf("electron with pt %4.1f\n", el.pt());
  selelecs.push_back(el);
  }
  */  
  return selelecs;
}

vector<cat::Jet> catJpsiAnalyzer::selectJets(const edm::View<cat::Jet>* jets, vector<TLorentzVector> recolep )
{
  vector<cat::Jet> seljets;
  for (auto jet : *jets) {
    if (!jet.LooseId()) continue;
    if (jet.pt() <= 30.) continue;
    if (fabs(jet.eta()) >= 2.4)	continue;
    if (jet.tlv().DeltaR(recolep[0]) <= 0.4) continue;
    if (jet.tlv().DeltaR(recolep[1]) <= 0.4) continue;
    // printf("jet with pt %4.1f\n", jet.pt());
    seljets.push_back(jet);
  }
  return seljets;
}

vector<cat::Jet> catJpsiAnalyzer::selectBJets(vector<cat::Jet> & jets )
{
  vector<cat::Jet> selBjets;
  if ( sysEnergy_ == 13 ) {
    for (auto jet : jets) {
      float jets_CSVInclV2 = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      if (jets_CSVInclV2 <= 0.814) continue;	
      //printf("b jet with pt %4.1f\n", jet.pt());
      selBjets.push_back(jet);
    }
  }
  else if (sysEnergy_ == 8 ) {
    for (auto jet : jets) {
      float jets_CSV = jet.bDiscriminator("combinedSecondaryVertexBJetTags");
      if ( jets_CSV <= 0.244) continue ;
      selBjets.push_back(jet); 
    }
  }
  else {
    std::cout<<"System Energy must be 8TeV or 13TeV at this script."<<std::endl;
  }
  return selBjets;
}

// ------------ method called once each job just before starting event loop  ------------
  void 
catJpsiAnalyzer::beginJob()
{
  //solver = new TtFullLepKinSolver(tmassbegin_, tmassend_, tmassstep_, nupars_);
}

// ------------ method called once each job just after ending the event loop  ------------
  void 
catJpsiAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
  void 
catJpsiAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
  void 
catJpsiAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
  void 
catJpsiAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
  void 
catJpsiAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
catJpsiAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(catJpsiAnalyzer);
