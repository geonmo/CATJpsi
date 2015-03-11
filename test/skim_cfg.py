import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIMMING")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
## Standard setup
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

## Options and Output Report
process.options = cms.untracked.PSet(
    #allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

## Source
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring())


## total event counter
process.totaEvents = cms.EDProducer("EventCountProducer")
process.p = cms.Path(process.totaEvents)

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('catTuple.root'),
    outputCommands = cms.untracked.vstring(
        'drop *_*_*_SKIMMING',
        'keep *',
    )
)
process.outpath = cms.EndPath(process.out)

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnMC: True default")
options.register('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "globalTag: 1  default")

options.parseArguments()
runOnMC = options.runOnMC
globalTag = options.globalTag

if globalTag:
    process.GlobalTag.globaltag = cms.string(globalTag)
if not globalTag:
    from Configuration.AlCa.autoCond import autoCond
    if runOnMC:
        process.GlobalTag.globaltag = autoCond['startup']
    else:
        process.GlobalTag.globaltag = autoCond['com10']
print "using globaltag", process.GlobalTag.globaltag

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('aodJpsiFiltered.root'),
    outputCommands = cms.untracked.vstring(
                                           'drop *_*_*_SKIMMING',
                                           'keep *',
    ), 
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('jpsifilter_path')
    )
)



## Max Number of Events
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
print options.maxEvents
#if ( options.maxEvents is not null ) : 
#  process.maxEvents.input = options.maxEvents

#process.source.fileNames = options.inputFiles
process.source.fileNames = cms.untracked.vstring('/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v2/10000/000560C1-FD97-E211-9F33-00304867924E.root')

## to suppress the long output at the end of the job
process.MessageLogger.cerr.threshold = ''
if options.maxEvents < 0:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = False

process.genJpsi = cms.EDFilter("PdgIdCandViewSelector",
    src = cms.InputTag("genParticles"),
    pdgId = cms.vint32( 443 )
)
process.jpsiFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("genJpsi"),
    minNumber = cms.uint32(1),
)

process.jpsifilter_path = cms.Path(process.genJpsi*process.jpsiFilter)
