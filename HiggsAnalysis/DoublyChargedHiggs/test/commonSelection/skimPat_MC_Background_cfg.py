import FWCore.ParameterSet.Config as cms

globalTag = 'MC_3XY_V26'

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('%s::All' % globalTag)
process.load("Configuration.StandardSequences.MagneticField_cff")

# Require HLT
process.load("HiggsAnalysis.DoublyChargedHiggs.hltFilters_cfi")

# HZZ Skimming
#process.load("HiggsAnalysis.Skimming.higgsToZZ4Leptons_Sequences_cff")
process.load("HiggsAnalysis.Skimming.higgsToZZ4Leptons_Filter_cfi")
process.higgsToZZ4Leptons_skimFilterSeq = cms.Sequence(process.higgsToZZ4LeptonsFilter)

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# Modification for Summer09-rereco sample
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
run33xOnReRecoMC( process, "ak5GenJets" )
#run33xOn31xMC(process)

# Good lepton selection
process.load("HiggsAnalysis.DoublyChargedHiggs.leptonSelectors_cff")

process.p = cms.Path(
    process.dhHLTFilters*
    process.higgsToZZ4Leptons_skimFilterSeq*
    process.patDefaultSequence*
    process.leptonSelectionSeq
)

# Output module configuration
from HiggsAnalysis.DoublyChargedHiggs.eventContent_cff import patDHEventContentWithLeptonFilter

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PAT.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    #save PAT Layer 1 output
    outputCommands = cms.untracked.vstring('drop *', *patDHEventContentWithLeptonFilter)
)
process.outpath = cms.EndPath(process.out)

