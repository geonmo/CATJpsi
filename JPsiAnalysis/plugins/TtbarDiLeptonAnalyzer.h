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

#include "CATJpsi/JPsiAnalysis/src/TMVAClassification_BDT.class.C"
using namespace std;

class TtbarDiLeptonAnalyzer : public edm::EDAnalyzer {
  public:
    explicit TtbarDiLeptonAnalyzer(const edm::ParameterSet&);
    ~TtbarDiLeptonAnalyzer();

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
    float passingSteps(int& trigger, int channel, float met, float ll_mass, float ll_charge, int selectedJets_size, int btag);
    const reco::Candidate* getLast(const reco::Candidate* p);

    TLorentzVector leafToTLorentzVector(reco::LeafCandidate & leaf)
    {return TLorentzVector(leaf.px(), leaf.py(),leaf.pz(),leaf.energy());}

    edm::InputTag muonToken_;
    edm::InputTag elecToken_;
    edm::InputTag jetToken_;
    edm::InputTag metToken_;
    edm::InputTag secToken_;
    edm::InputTag vtxToken_;
    edm::InputTag mcLabel_;
    edm::InputTag partonTop_channel_;
    edm::InputTag partonTop_modes_;

    int sysEnergy_;

    TTree * ttree_;
    int b_genChannel, b_genMode1, b_genMode2, b_partonChannel, b_partonMode1, b_partonMode2;
    int b_njet, b_nbjet, b_step, b_channel, b_inPhase;
    float b_MET, b_maxweight;

    float b_lep1_pt, b_lep1_eta, b_lep1_phi;
    float b_lep2_pt, b_lep2_eta, b_lep2_phi;
    float b_ll_pt, b_ll_eta, b_ll_phi, b_ll_m;

    int b_pvN;
    float b_pileupWeight, b_pileupWeight_up, b_pileupWeight_dn;      


    typedef std::vector<float> floats;
    typedef std::vector<int> ints;
    floats * b_jpsi_pt, * b_jpsi_eta, * b_jpsi_phi, * b_jpsi_mass;
    floats * b_jpsi_vProb, * b_jpsi_lxy, * b_jpsi_dca, * b_jpsi_minDR, * b_jpsi_minBDR, * b_jpsi_mva;
    ints * b_jpsi_muID, * b_jpsi_trackQuality;

    floats * b_ljpsi1_pt, * b_ljpsi1_eta, * b_ljpsi1_phi, * b_ljpsi1_m;
    floats * b_ljpsi2_pt, * b_ljpsi2_eta, * b_ljpsi2_phi, * b_ljpsi2_m;

    TtFullLepKinSolver* solver;
    double tmassbegin_, tmassend_, tmassstep_;
    vector<double> nupars_;

    bool runOnMC_;
    int gen_channel;
    ReadBDT* bdt_; 
    std::vector<int> gen_modes;
    //enum TTbarMode { CH_NONE = 0, CH_FULLHADRON = 1, CH_SEMILEPTON, CH_FULLLEPTON };
    //enum DecayMode { CH_HADRON = 1, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };
};
