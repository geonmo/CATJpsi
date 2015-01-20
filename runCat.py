from PhysicsTools.PatAlgos.patTemplate_cfg import *
#from CATTools.CatProducer.catTemplate_cfg import *
## some options
doSecVertex=True # for jpsi candidates
    
## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnMC: True default")

options.parseArguments()
runOnMC = options.runOnMC

####################################################################################################
## running PAT
postfix = "PFlow"
jetAlgo="AK5"
from CATTools.CatProducer.catPatSetup_cff import *
catPatConfig(process, runOnMC, postfix, jetAlgo)

process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Winter14_V5_DATA_AK5PF'),
            ),
      ), 
      connect = cms.string('sqlite:Winter14_V5_DATA.db')
)

process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')
####################################################################################################

from CATTools.CatProducer.catSetup_cff import *
catSetup(process, runOnMC, doSecVertex)

process.maxEvents.input = options.maxEvents

process.source.fileNames = options.inputFiles

## to suppress the long output at the end of the job
process.MessageLogger.cerr.threshold = ''
if options.maxEvents < 0:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = False

# Remove pt>5GeV/c cut for Muon.
process.pfSelectedMuonsPFlow.cut = cms.string('')
# Jet Correction option for 2012Rereco 
process.pfPileUpPFlow.checkClosestZVertex = False

# top projections in PF2PAT: we are turning off top projection
process.pfNoMuonPFlow.enable = False #True
process.pfNoElectronPFlow.enable = False #True
process.pfNoJetPFlow.enable = False #True

# Use non-isolated muons and electrons
process.patMuonsPFlow.pfMuonSource = "pfMuonsPFlow"
process.patElectronsPFlow.pfElectronSource = "pfElectronsPFlow"
# And turn on delta-beta corrections while building pfIsolated*PFlow
process.pfIsolatedMuonsPFlow.doDeltaBetaCorrection = True
process.pfIsolatedElectronsPFlow.doDeltaBetaCorrection = True


