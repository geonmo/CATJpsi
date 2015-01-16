from PhysicsTools.PatAlgos.patTemplate_cfg import *
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register ('runOnMC', 1,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "runOnMC")
import sys
if hasattr(sys, "argv") == True:
    options.parseArguments()
    runOnMC = options.runOnMC

postfix = "PFlow"
jetAlgo="AK5"
doSecVertex=True # for jpsi candidates

from CATTools.CatProducer.catPatSetup_cff import *
from CATTools.CatProducer.catSetup_cff import *
catPatConfig(process, runOnMC, postfix, jetAlgo)
catSetup(process, runOnMC, doSecVertex)

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
process.maxEvents.input = -1
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source.fileNames = cms.untracked.vstring(
'file:$HOME/jpsi.root'
)


addGenParticle = cms.untracked.vstring('keep recoGenParticles_genParticles*_*_*')
process.out.outputCommands += addGenParticle

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


