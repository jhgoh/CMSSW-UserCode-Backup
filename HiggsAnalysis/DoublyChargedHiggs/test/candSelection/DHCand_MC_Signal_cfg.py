import FWCore.ParameterSet.Config as cms

# Set variables from the os environment
globalTag = 'MC_3XY_V26'

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
process.load("HiggsAnalysis.DoublyChargedHiggs.genEventFilters_cff")
process.load("HiggsAnalysis.DoublyChargedHiggs.leptonSelectors_cff")
process.load("HiggsAnalysis.DoublyChargedHiggs.compositeCandProducers_cff")

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('DHCand.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('keep *')
)
process.outpath = cms.EndPath(process.out)

## Set run path
#analysisMode = 'MMMM'
#analysisMode = 'EEEE'
analysisMode = 'EMEM'

process.dhEventFilter_step = cms.Sequence(getattr(process, 'dhGenEventFilterTo'+analysisMode))

leptonSelectors = None
if 'M' in analysisMode:
    if leptonSelectors == None:
        leptonSelectors = process.goodPatMuons
    else:
        leptonSelectors += process.goodPatMuons
if 'E' in analysisMode:
    if leptonSelectors == None:
        leptonSelectors = process.goodPatElectrons
    else:
        leptonSelectors += process.goodPatElectrons

process.goodLeptonSelector_step = cms.Sequence(leptonSelectors)

candProducerModules = None
for i in set( (analysisMode[0:2], analysisMode[2:4]) ):
    if candProducerModules == None:
        candProducerModules = getattr(process, 'dhCandProducerTo'+i)
    else:
        candProducerModules += getattr(process, 'dhCandProducerTo'+i)
process.dhCandProducer_step = cms.Sequence(candProducerModules)

process.p = cms.Path(process.dhEventFilter_step*process.dhCandProducer_step)
