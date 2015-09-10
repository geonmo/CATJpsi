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
#include "CATTools/DataFormats/interface/SecVertex.h"

#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

class catJpsiAnalyzer : public edm::EDAnalyzer {
  public:
    explicit catJpsiAnalyzer(const edm::ParameterSet&);
    ~catJpsiAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    vector<cat::Muon> selectMuons(const edm::View<cat::Muon>* muons );
    vector<cat::Electron> selectElecs(const edm::View<cat::Electron>* elecs );
    vector<cat::Jet> selectJets(const edm::View<cat::Jet>* jets, vector<TLorentzVector> recolep);
    vector<cat::Jet> selectBJets(vector<cat::Jet> & jets );
    const reco::Candidate* getLast(const reco::Candidate* p);
    bool isEqual( TLorentzVector& gen, TLorentzVector& reco, bool exact);

    TLorentzVector leafToTLorentzVector(reco::LeafCandidate & leaf)
    {return TLorentzVector(leaf.px(), leaf.py(),leaf.pz(),leaf.energy());}

    edm::InputTag muonToken_;
    edm::InputTag elecToken_;
    edm::InputTag jetToken_;
    edm::InputTag metToken_;
    edm::InputTag secToken_;
    edm::InputTag vtxToken_;
    edm::InputTag mcLabel_;

    int sysEnergy_;

    TTree * ttree_;
    TTree * matchedtree_;
    TTree * unmatchedtree_;

    typedef std::vector<float> floats;
    typedef std::vector<int> ints;
    float b_lep1_pt, b_lep1_eta, b_lep1_phi, b_lep1_m;
    float b_lep2_pt, b_lep2_eta, b_lep2_phi, b_lep2_m;

    float  b_mjpsi_pt,     b_mjpsi_eta,  b_mjpsi_phi,  b_mjpsi_m;
    float  b_mjpsi_vProb,  b_mjpsi_lxy,  b_mjpsi_dca,  b_mjpsi_minJetDR, b_mjpsi_minBJetDR;
    int    b_mjpsi_muID,   b_mjpsi_trackQuality;
    float  b_mjpsi_minGenDR;

    float  b_ujpsi_pt,     b_ujpsi_eta,  b_ujpsi_phi,  b_ujpsi_m;
    float  b_ujpsi_vProb,  b_ujpsi_lxy,  b_ujpsi_dca,  b_ujpsi_minJetDR, b_ujpsi_minBJetDR;
    int    b_ujpsi_muID,   b_ujpsi_trackQuality;
    TtFullLepKinSolver* solver;

    bool runOnMC_;
};
