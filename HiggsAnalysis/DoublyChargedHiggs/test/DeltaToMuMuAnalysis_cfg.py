import FWCore.ParameterSet.Config as cms

# Set variables from the os environment
import os
CMSSWVersion = os.environ['CMSSW_VERSION']
dataTier = 'GEN-SIM-DIGI-RECO'
#globalTag = 'STARTUP_V12'
globalTag = 'IDEAL_V12'
if 'SAMPLENAME' in os.environ:
        sampleName = os.environ['SAMPLENAME']
else:
	sampleName = 'Hpp130MuMu_FastSim'

# Load Standard CMSSW process initial configurations
process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = cms.string('%s::All' % globalTag)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
)

# Set datafiles
process.source.fileNames.append('file:/pnfs/jhgoh/DoublyChargedHiggs/%s/%s_%s.root' % (CMSSWVersion, sampleName, globalTag))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# PAT sequences
process.load("PhysicsTools/PatAlgos/patSequences_cff")
process.load("PhysicsTools/PatAlgos/patEventContent_cff")
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

# Muon selection
process.posMuons = cms.EDProducer("PATMuonCandSelector",
    muon = cms.InputTag("cleanLayer1Muons"),
    muonCutSet = cms.PSet(
        charge = cms.untracked.int32(1),
        useGlobalMuonsOnly = cms.untracked.bool(True),
        minPt = cms.untracked.double(7)
    )
)

process.negMuons = cms.EDProducer("PATMuonCandSelector",
    muon = cms.InputTag("cleanLayer1Muons"),
    muonCutSet = cms.PSet(
        charge = cms.untracked.int32(-1),
        useGlobalMuonsOnly = cms.untracked.bool(True),
        minPt = cms.untracked.double(7)
    )
)

# Candidate production
process.posDeltaToMuMu = cms.EDProducer("DimuonProducer",
    muon1 = cms.InputTag("posMuons"),
    muon2 = cms.InputTag("posMuons"),
    fitterType = cms.string("None"), #KalmanVertexFitter"),
    vertexFitSet = cms.PSet(
        maxDistance = cms.double(10.0)
    )
)

process.negDeltaToMuMu = cms.EDProducer("DimuonProducer",
    muon1 = cms.InputTag("negMuons"),
    muon2 = cms.InputTag("negMuons"),
    fitterType = cms.string("iNone"), #KalmanVertexFitter"),
    vertexFitSet = cms.PSet(
        maxDistance = cms.double(10.0)
    )
)

# User analysis block
process.fourMuonAnalyzer = cms.EDAnalyzer("FourMuonAnalyzer",
    posDelta = cms.InputTag("posDeltaToMuMu"),
    negDelta = cms.InputTag("negDeltaToMuMu"),
    deltaCutSet = cms.PSet(
        maxNormalizedChi2 = cms.untracked.double(1000),
        minPt = cms.untracked.double(7)
    ),
    nInterested = cms.untracked.uint32(2)
)

# File output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('h_%s_%s.root' % (sampleName, globalTag))
)

# Module sequences and Paths
process.posDeltaSeq = cms.Sequence(process.posMuons*process.posDeltaToMuMu)
process.negDeltaSeq = cms.Sequence(process.negMuons*process.negDeltaToMuMu)
process.deltaCombineSeq = cms.Sequence(process.posDeltaSeq+process.negDeltaSeq)
process.deltaAnalysisSeq = cms.Sequence(process.fourMuonAnalyzer)

process.p = cms.Path(process.deltaCombineSeq*
                     process.deltaAnalysisSeq)
