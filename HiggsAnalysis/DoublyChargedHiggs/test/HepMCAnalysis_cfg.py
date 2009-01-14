import FWCore.ParameterSet.Config as cms
            
process = cms.Process('USER')

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
	fileNames = cms.untracked.vstring('file:Hpp160MuMu_cfi_py_GEN.root')
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.hepMCAnalyzer = cms.EDAnalyzer("Hpp2MuHepMCAnalyzer")
process.analyzer_step = cms.Path(process.hepMCAnalyzer)
process.schedule = cms.Schedule(process.analyzer_step);
