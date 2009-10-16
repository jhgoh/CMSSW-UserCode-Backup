import FWCore.ParameterSet.Config as cms

globalTag = 'MC_31X_V8'

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('PATSummaryTables')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
#    PATSummaryTables = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
#)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

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

# Set datafiles
process.source.fileNames.append('file:res/step1.root')
process.source.secondaryFileNames.append('file:res/step2.root')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.p = cms.Path(
    process.patDefaultSequence  
)

# Output module configuration
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
patEventContent = [
    'keep *_cleanLayer1Photons_*_*', 
    'keep *_cleanLayer1Electrons_*_*', 
    'keep *_cleanLayer1Muons_*_*', 
    'keep *_cleanLayer1Taus_*_*', 
    'keep *_cleanLayer1Jets_*_*',
    'keep *_layer1METs_*_*',
    'keep *_cleanLayer1Hemispheres_*_*',
    'keep *_cleanLayer1PFParticles_*_*',
    'keep *_offlinePrimaryVertices_*_*',
    'keep recoGsfTrackExtras_pixelMatchGsfFit_*_*',
    'keep recoTrackExtras_pixelMatchGsfFit_*_*',
    'keep recoTracks_generalTracks_*_*',
    'keep recoTrackExtras_generalTracks_*_*',
    'keep *_offlineBeamSpot_*_*'
]

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PAT.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output
    outputCommands = cms.untracked.vstring('drop *', *patEventContent ) # you need a '*' to unpack the list of commands 'patEventContent'
    #outputCommands = cms.untracked.vstring('keep *')
)
process.outpath = cms.EndPath(process.out)

