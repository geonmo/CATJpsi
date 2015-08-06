import FWCore.ParameterSet.Config as cms

process = cms.Process("TtbarDiLeptonAnalyzer")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.SkipEvent = cms.untracked.vstring('ProductNotFound')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.partonTop = cms.EDProducer("PartonTopProducer",genParticles = cms.InputTag("prunedGenParticles"))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(),
)
#for i in xrange(1,101):
#    process.source.fileNames.append('file:/cms/data/xrd/store/user/jlee/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/cat74v2_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150713_164609/0000/catTuple_%d.root' % i)

#process.source.fileNames.append('file:/pnfs/user/jlee/scratch/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/catTuple.root')
#process.source.fileNames.append('file:/store/group/CAT/SingleMuon/v7-3-0_Run2015B-PromptReco-v1/150720_060727/0000/catTuple_1.root')
process.source.fileNames.append('/store/user/geonmo/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola//catJpsi_20150805_Summer12_DR53X-PU_S10_START53_V7C-v2/150805_155516/0000/catTuple_983.root')

process.ttll = cms.EDAnalyzer("TtbarDiLeptonAnalyzer",
    vertices = cms.InputTag("goodOfflinePrimaryVertices"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    mcLabel = cms.InputTag("genParticles"),
    partonTop_channel = cms.InputTag("partonTop","channel"),
    partonTop_modes = cms.InputTag("partonTop", "modes"),

    tmassbegin = cms.double(100),
    tmassend   = cms.double(300),
    tmassstep  = cms.double(  1),
    neutrino_parameters = cms.vdouble(27.23,53.88,19.92,53.89,19.9)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("top.root"
))

process.p = cms.Path(process.ttll)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
