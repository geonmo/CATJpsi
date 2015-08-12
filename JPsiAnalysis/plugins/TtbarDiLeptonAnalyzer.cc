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
#include "CATJpsi/JPsiAnalysis/plugins/TtbarDiLeptonAnalyzer.h"

#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

TtbarDiLeptonAnalyzer::TtbarDiLeptonAnalyzer(const edm::ParameterSet& iConfig)
{
  muonToken_ = iConfig.getParameter<edm::InputTag>("muons");
  elecToken_ = iConfig.getParameter<edm::InputTag>("electrons");
  jetToken_  = iConfig.getParameter<edm::InputTag>("jets");
  metToken_  = iConfig.getParameter<edm::InputTag>("mets");     
  secToken_  = iConfig.getParameter<edm::InputTag>("secVtxs");     
  vtxToken_  = iConfig.getParameter<edm::InputTag>("vertices");
  mcLabel_   = iConfig.getParameter<edm::InputTag>("mcLabel");
  partonTop_channel_ = iConfig.getParameter<edm::InputTag>("partonTop_channel");
  partonTop_modes_   = iConfig.getParameter<edm::InputTag>("partonTop_modes");

  tmassbegin_     = iConfig.getParameter<double>       ("tmassbegin");
  tmassend_       = iConfig.getParameter<double>       ("tmassend");
  tmassstep_      = iConfig.getParameter<double>       ("tmassstep");
  nupars_         = iConfig.getParameter<vector<double> >("neutrino_parameters");

  edm::Service<TFileService> fs;
  ttree_ = fs->make<TTree>("tree", "tree");
  ttree_->Branch("gen_channel", &b_genChannel, "gen_channel/I");
  ttree_->Branch("gen_mode1", &b_genMode1, "gen_mode1/I");
  ttree_->Branch("gen_mode2", &b_genMode2, "gen_mode2/I");
  ttree_->Branch("parton_channel", &b_partonChannel, "parton_channel/I");
  ttree_->Branch("parton_mode1", &b_partonMode1, "parton_mode1/I");
  ttree_->Branch("parton_mode2", &b_partonMode2, "parton_mode2/I");

  ttree_->Branch("njet", &b_njet, "njet/I");
  ttree_->Branch("nbjet", &b_nbjet, "nbjet/I");
  ttree_->Branch("MET", &b_MET, "MET/F");
  ttree_->Branch("channel", &b_channel, "channel/I");
  ttree_->Branch("step", &b_step, "step/I");
  ttree_->Branch("inPhase", &b_inPhase, "inPhase/I");

  ttree_->Branch("lep1_pt", &b_lep1_pt, "lep1_pt/F");
  ttree_->Branch("lep1_eta", &b_lep1_eta, "lep1_eta/F");
  ttree_->Branch("lep1_phi", &b_lep1_phi, "lep1_phi/F");
  ttree_->Branch("lep2_pt", &b_lep2_pt, "lep2_pt/F");
  ttree_->Branch("lep2_eta", &b_lep2_eta, "lep2_eta/F");
  ttree_->Branch("lep2_phi", &b_lep2_phi, "lep2_phi/F");
  ttree_->Branch("ll_pt", &b_ll_pt, "ll_pt/F");
  ttree_->Branch("ll_eta", &b_ll_eta, "ll_eta/F");
  ttree_->Branch("ll_phi", &b_ll_phi, "ll_phi/F");
  ttree_->Branch("ll_m", &b_ll_m, "ll_m/F");

  b_jpsi_pt = new floats;
  b_jpsi_eta = new floats;
  b_jpsi_phi = new floats;
  b_jpsi_mass = new floats;
  b_jpsi_vProb = new floats;
  b_jpsi_l3D = new floats;
  b_jpsi_dca = new floats;
  b_jpsi_muID = new floats;
  b_jpsi_trackQuality = new floats;
  b_jpsi_minDR = new floats;
  b_jpsi_minBDR = new floats;

  ttree_->Branch("jpsi_pt",&b_jpsi_pt);
  ttree_->Branch("jpsi_eta",&b_jpsi_eta);
  ttree_->Branch("jpsi_phi",&b_jpsi_phi);
  ttree_->Branch("jpsi_mass",&b_jpsi_mass);
  ttree_->Branch("jpsi_vProb",&b_jpsi_vProb);
  ttree_->Branch("jpsi_l3D",&b_jpsi_l3D);
  ttree_->Branch("jpsi_dca",&b_jpsi_dca);
  ttree_->Branch("jpsi_muID",&b_jpsi_muID);
  ttree_->Branch("jpsi_trackQuality",&b_jpsi_trackQuality);
  ttree_->Branch("jpsi_minDR",&b_jpsi_minDR);
  ttree_->Branch("jpsi_minBDR",&b_jpsi_minBDR);

}
TtbarDiLeptonAnalyzer::~TtbarDiLeptonAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void TtbarDiLeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  b_genChannel = -1; 
  b_genMode1 = -1;
  b_genMode2 = -1; 
  b_partonChannel = -1; 
  b_partonMode1 = -1; 
  b_partonMode2 = -1; 
  b_MET = -1; 
  b_njet = -1;
  b_nbjet = -1;
  b_channel = -1;
  b_step = -1;
  b_inPhase = 0;
  b_lep1_pt = -9; b_lep1_eta = -9; b_lep1_phi = -9;
  b_lep2_pt = -9; b_lep2_eta = -9; b_lep2_phi = -9;
  b_ll_pt = -9; b_ll_eta = -9; b_ll_phi = -9; b_ll_m = -9;
  b_jpsi_pt->clear();
  b_jpsi_eta->clear();
  b_jpsi_phi->clear();
  b_jpsi_mass->clear();
  b_jpsi_vProb->clear();
  b_jpsi_l3D->clear();
  b_jpsi_dca->clear();
  b_jpsi_muID->clear();
  b_jpsi_trackQuality->clear();
  b_jpsi_minDR->clear();
  b_jpsi_minBDR->clear();
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
  //edm::Handle<int> partonTop_channel;
  //edm::Handle<vector<int> > partonTop_modes;

  edm::Handle<edm::View<cat::SecVertex> > svertxs ;
  iEvent.getByLabel(secToken_, svertxs);

  if (runOnMC_){
    int nMuon = 0;
    int nElectron = 0;
    gen_modes.clear();

    iEvent.getByLabel(mcLabel_,genParticles); 
    for (const reco::GenParticle & g : *genParticles){
      const reco::Candidate* w=0;
      const reco::Candidate* wLast=0;    
      const reco::Candidate* lep=0;
      if (fabs(g.pdgId()) == 6){ 
        for (unsigned int i = 0; i < g.numberOfDaughters(); ++i){
          if (fabs(g.daughter(i)->pdgId())  == 24){ w = g.daughter(i); break; }
        }
      }
      if (w){
        wLast=getLast(w);
        for (unsigned int i = 0; i < wLast->numberOfDaughters(); ++i){
          if ((fabs(wLast->daughter(i)->pdgId()) == 11) || (fabs(wLast->daughter(i)->pdgId()) == 13) || (fabs(wLast->daughter(i)->pdgId()) == 15)){
            lep = wLast->daughter(i);
            break;
          }
        }
      }

      if (lep){
        int mode = 1;
        if ( fabs(lep->pdgId()) == 13){ ++nMuon; mode = 2; }
        else if ( fabs(lep->pdgId()) == 11){ ++nElectron; mode = 3; }
        else if ( fabs(lep->pdgId()) == 15){
          for (unsigned int i = 0; i < lep->numberOfDaughters(); ++i){
            if ( fabs(lep->daughter(i)->pdgId()) == 13 ) { mode = 5; break;}
            else if ( fabs(lep->daughter(i)->pdgId()) == 11 ) { mode = 6; break;}
            mode = 4;
          }
        }
        gen_modes.push_back(mode);
      }
    }

    if ( gen_modes.size() == 0 ) { gen_modes.push_back(0); }
    if ( gen_modes.size() == 1 ) { gen_modes.push_back(0); }

    gen_channel = 0;
    const int nLepton = nElectron + nMuon;
    gen_channel = nLepton+1;

    b_genChannel = gen_channel; 
    b_genMode1 = gen_modes[0]; 
    b_genMode2 = gen_modes[1]; 

  }
  vector<cat::Muon> selectedMuons = selectMuons( muons.product() );
  vector<cat::Electron> selectedElectrons = selectElecs( electrons.product() );

  vector<TLorentzVector> recolep; 
  for (auto lep : selectedMuons){ recolep.push_back(lep.tlv()); }
  for (auto lep : selectedElectrons){ recolep.push_back(lep.tlv()); }
  if (recolep.size() != 2){
    ttree_->Fill();
    return;
  }

  b_lep1_pt = recolep[0].Pt();
  b_lep1_eta = recolep[0].Eta();
  b_lep1_phi = recolep[0].Phi();
  b_lep2_pt = recolep[1].Pt();
  b_lep2_eta = recolep[1].Eta();
  b_lep2_phi = recolep[1].Phi();

  if ((recolep[0].Pt() > 20) && (recolep[1].Pt() > 20) && (fabs(recolep[0].Eta()) < 2.4) && (fabs(recolep[1].Eta()) < 2.4)) b_inPhase = 1;

  float channel = selectedElectrons.size();

  float ll_charge = 0. ;
  if (channel == 0) ll_charge = selectedMuons[0].charge()*selectedMuons[1].charge();
  if (channel == 1) ll_charge = selectedMuons[0].charge()*selectedElectrons[0].charge();
  if (channel == 2) ll_charge = selectedElectrons[0].charge()*selectedElectrons[1].charge();

  vector<cat::Jet> selectedJets = selectJets( jets.product(), recolep );
  vector<cat::Jet> selectedBJets = selectBJets( selectedJets );

  //  printf("selectedMuons %lu, selectedElectrons %lu, selectedJets %lu, selectedBJets %lu\n",selectedMuons.size(), selectedElectrons.size(), selectedJets.size(), selectedBJets.size() );

  TLorentzVector met = mets->front().tlv();

  TLorentzVector tlv_ll = recolep[0]+recolep[1];
  b_ll_pt = tlv_ll.Pt();
  b_ll_eta = tlv_ll.Eta();
  b_ll_phi = tlv_ll.Phi();
  b_ll_m = tlv_ll.M();

  b_MET = met.Pt(); 
  b_njet = selectedJets.size();
  b_nbjet = selectedBJets.size();
  b_channel = channel;

  auto svtxs = svertxs.product();
  for ( auto catJpsi = svtxs->begin(), end= svtxs->end() ; catJpsi != end ; ++catJpsi) {
    b_jpsi_pt->push_back(catJpsi->pt());
    b_jpsi_eta->push_back(catJpsi->eta());
    b_jpsi_phi->push_back(catJpsi->phi());
    b_jpsi_mass->push_back(catJpsi->mass());
    b_jpsi_vProb->push_back(catJpsi->vProb());
    b_jpsi_l3D->push_back(catJpsi->l3D());
    b_jpsi_dca->push_back(catJpsi->dca());
    b_jpsi_muID->push_back(catJpsi->muID());
    b_jpsi_trackQuality->push_back(catJpsi->trackQuality());
    float min_DR = 999;
    float min_BDR = 999;

    for (auto jet = selectedJets.begin(), end = selectedJets.end(); jet != end; ++jet){
      float csv = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
      TLorentzVector jpsi_tlv = TLorentzVector(catJpsi->px(), catJpsi->py(),catJpsi->pz(),catJpsi->energy()); 
      TLorentzVector jet_tlv = jet->tlv();
      float deltaR = jpsi_tlv.DeltaR( jet_tlv) ;
      if ( min_DR > deltaR ) min_DR = deltaR; 
      if ( csv > 0.244 && min_BDR > deltaR) min_BDR = deltaR; 
    }
    b_jpsi_minDR->push_back( min_DR ); 
    b_jpsi_minBDR->push_back( min_BDR ); 
  }
  float step = passingSteps( channel, met.Pt(), (recolep[0]+recolep[1]).M(), ll_charge, selectedJets.size(), selectedBJets.size() );
  b_step = step;

  ////////////////////////////////////////////////////////  KIN  /////////////////////////////////////
  int kin=0; TLorentzVector nu1, nu2, top1, top2;
  double maxweight=0;
  cat::Jet kinj1, kinj2;

  for (auto jet1 = selectedJets.begin(), end = selectedJets.end(); jet1 != end; ++jet1){
    for (auto jet2 = next(jet1); jet2 != end; ++jet2){

      double weight1 =0; double weight2 =0;
      TLorentzVector nu11, nu12, nu21, nu22;
      TLorentzVector recojet1= jet1->tlv();
      TLorentzVector recojet2= jet2->tlv();

      double xconstraint = recolep[0].Px()+recolep[1].Px()+ (recojet1).Px() + (recojet2).Px() +met.Px();
      double yconstraint = recolep[0].Py()+recolep[1].Py()+ (recojet2).Py() + (recojet1).Py() +met.Py();

      solver->SetConstraints(xconstraint, yconstraint);
      TtFullLepKinSolver::NeutrinoSolution nuSol= solver->getNuSolution( recolep[0], recolep[1] , recojet1, recojet2);
      weight1 = nuSol.weight;
      nu11 = leafToTLorentzVector(nuSol.neutrino);
      nu12 = leafToTLorentzVector(nuSol.neutrinoBar);

      TtFullLepKinSolver::NeutrinoSolution nuSol2= solver->getNuSolution( recolep[0], recolep[1] , recojet2, recojet1);
      weight2 = nuSol2.weight;
      nu21 = leafToTLorentzVector(nuSol2.neutrino);
      nu22 = leafToTLorentzVector(nuSol2.neutrinoBar);
      if (weight1 > maxweight || weight2 > maxweight){
        if(weight1>weight2 && weight1>0){
          maxweight = weight1; kinj1=(*jet1); kinj2=(*jet2); nu1 = nu11; nu2 = nu12; kin++;
          top1 = recolep[0]+recojet1+nu11; top2 = recolep[1]+recojet2+nu12;
        }
        else if(weight2>weight1 && weight2>0){
          maxweight = weight2; kinj1=(*jet2); kinj2=(*jet1); nu1 = nu21; nu2 = nu22; kin++;
          top1 = recolep[0]+recojet2+nu21; top2 = recolep[1]+recojet1+nu22;
        }
      }
    }
  }
  b_maxweight = maxweight;
  //  printf("maxweight %f, top1.M() %f, top2.M() %f \n",maxweight, top1.M(), top2.M() );
  // printf("%2d, %2d, %2d, %2d, %6.2f, %6.2f, %6.2f\n", b_njet, b_nbjet, b_step, b_channel, b_MET, b_ll_mass, b_maxweight);

  ttree_->Fill();
}

const reco::Candidate* TtbarDiLeptonAnalyzer::getLast(const reco::Candidate* p)
{
  for ( size_t i=0, n=p->numberOfDaughters(); i<n; ++i )
  {
    const reco::Candidate* dau = p->daughter(i);
    if ( p->pdgId() == dau->pdgId() ) return getLast(dau);
  }
  return p;
}

vector<cat::Muon> TtbarDiLeptonAnalyzer::selectMuons(const edm::View<cat::Muon>* muons )
{
  vector<cat::Muon> selmuons;
  for (auto mu : *muons) {
    //if (!mu.isMediumMuon()) continue;
    if (!mu.isTightMuon()) continue;
    if (mu.pt() <= 20.) continue;
    if (fabs(mu.eta()) >= 2.4) continue;
    if (mu.relIso(0.4) >= 0.12) continue;
    //printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    selmuons.push_back(mu);
  }

  return selmuons;
}

vector<cat::Electron> TtbarDiLeptonAnalyzer::selectElecs(const edm::View<cat::Electron>* elecs )
{
  vector<cat::Electron> selelecs;
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
  return selelecs;
}

vector<cat::Jet> TtbarDiLeptonAnalyzer::selectJets(const edm::View<cat::Jet>* jets, vector<TLorentzVector> recolep )
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

vector<cat::Jet> TtbarDiLeptonAnalyzer::selectBJets(vector<cat::Jet> & jets )
{
  vector<cat::Jet> selBjets;
  for (auto jet : jets) {
    float jets_CSVInclV2 = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    if (jets_CSVInclV2 <= 0.814) continue;	
    //printf("b jet with pt %4.1f\n", jet.pt());
    selBjets.push_back(jet);
  }
  return selBjets;
}

float TtbarDiLeptonAnalyzer::passingSteps(int channel, float met, float ll_mass, float ll_charge, int selectedJets_size, int btag)
{
  int step = 0;
  if (ll_mass <= 20.) return step;
  if (ll_charge > 0.) return step;
  step = 1;
  if (channel != 1){
    if ((ll_mass > 76) and (ll_mass < 106)) return step;
  }
  step = 2;
  if (selectedJets_size < 2) return step;
  step = 3;
  if (channel == 1){
    step = 4;
  }
  else{
    if (met <= 40.) return step;
  }
  step = 4;
  if (btag <= 0) return step;
  step = 5;

  return step;
}

// ------------ method called once each job just before starting event loop  ------------
  void 
TtbarDiLeptonAnalyzer::beginJob()
{
  solver = new TtFullLepKinSolver(tmassbegin_, tmassend_, tmassstep_, nupars_);
}

// ------------ method called once each job just after ending the event loop  ------------
  void 
TtbarDiLeptonAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
  void 
TtbarDiLeptonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
  void 
TtbarDiLeptonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
  void 
TtbarDiLeptonAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
  void 
TtbarDiLeptonAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TtbarDiLeptonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TtbarDiLeptonAnalyzer);
