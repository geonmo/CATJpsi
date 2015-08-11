import FWCore.ParameterSet.Config as cms

process = cms.Process("Matching")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.source.fileNames = cms.untracked.vstring('/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v2/10001/EEECECEE-DF98-E211-B5D0-0025905964A2.root')
#process.source.fileNames = options.inputFiles
#process.source.skipEvents=cms.untracked.uint32(7000)

## to suppress the long output at the end of the job

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(),)

process.source.fileNames.append('/store/user/geonmo/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola//catJpsi_20150805_Summer12_DR53X-PU_S10_START53_V7C-v2/150805_155516/0000/catTuple_983.root')

process.MessageLogger.cerr.threshold = ''


process.SecVertexsGenParticleMatch = cms.EDProducer("MCMatcher",     # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src     = cms.InputTag("catSecVertexs"),        # RECO objects to match
    matched = cms.InputTag("genParticles"), # mc-truth particle collection
    mcPdgId     = cms.vint32(443),           # one or more PDG ID (13 = muon); absolute values (see below)
    checkCharge = cms.bool(True),           # True = require RECO and MC objects to have the same charge
    mcStatus = cms.vint32(1),               # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR = cms.double(0.15),            # Minimum deltaR for the match
    maxDPtRel = cms.double(0.05),            # Minimum deltaPt/Pt for the match
    resolveAmbiguities = cms.bool(True),    # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True) # False = just match input in order; True = pick lowest deltaR pair first
)

process.p = cms.Path(process.SecVertexsGenParticleMatch)
process.out =  cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('store.root'),
    outputCommands = cms.untracked.vstring('keep *',),
)
process.outpath = cms.EndPath(process.out)
