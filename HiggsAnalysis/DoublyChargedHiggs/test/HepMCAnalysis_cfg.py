import FWCore.ParameterSet.Config as cms
import os
            
process = cms.Process('USER')
higgsMass = os.environ['HIGGSMASS']

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedEarly10TeVCollision_cff')           
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.GlobalTag.globaltag = "IDEAL_V9::All"

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring('file:Hpp%sMuMu_cfi_py_GEN.root' % higgsMass)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('Result_%s.root' % higgsMass)
)

process.hepMCAnalyzer = cms.EDAnalyzer("Hpp2MuHepMCAnalyzer",
    nBinPt = cms.untracked.uint32(50),
    nBinEta = cms.untracked.uint32(100),
    nBinM = cms.untracked.uint32(50),

    minTrkPt = cms.untracked.double(0),
    maxTrkPt = cms.untracked.double(500),

    minHiggsPt = cms.untracked.double(0),
    maxHiggsPt = cms.untracked.double(100),

    minEta = cms.untracked.double(-2.5),
    maxEta = cms.untracked.double(2.5),

    minM = cms.untracked.double(120),
    maxM = cms.untracked.double(320)
)
process.analyzer_step = cms.Path(process.hepMCAnalyzer)
process.schedule = cms.Schedule(process.analyzer_step);

# vim:set ts=4 sts=4 sw=4 expandtab: #
