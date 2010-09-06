import FWCore.ParameterSet.Config as cms

# Set variables from the os environment
globalTag = 'START3X_V26'

# Load Standard CMSSW process initial configurations
process = cms.Process("DHCand")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = cms.string('%s::All' % globalTag)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
)

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Load modules
process.load("HiggsAnalysis.DoublyChargedHiggs.compositeCandProducers_cff")

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('DHCand.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('keep *')
)
process.outpath = cms.EndPath(process.out)

## Set decay mode
#analysisMode = 'MMMM'
#analysisMode = 'EEEE'
analysisMode = 'EMEM'

candProducerModules = None
for i in set( (analysisMode[0:2], analysisMode[2:4]) ):
    if candProducerModules == None:
        candProducerModules = getattr(process, 'dhCandTo'+i)
    else:
        candProducerModules += getattr(process, 'dhCandTo'+i)
process.dhCand_step = cms.Sequence(candProducerModules)

process.p = cms.Path(
    process.dhCand_step * 
    process.atLeastOneDHCandFilter
)
