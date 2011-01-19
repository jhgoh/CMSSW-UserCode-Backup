# This config file is prepared to produce new streamer file (*.dat). 
# In outputCommands, you can put whatever hlt products you want to put in.

import FWCore.ParameterSet.Config as cms

process = cms.Process("TRANSFER")

import FWCore.Framework.test.cmsExceptionsFatal_cff
process.options = FWCore.Framework.test.cmsExceptionsFatal_cff.options

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:inputfile.root'
    )
)

process.a1 = cms.EDAnalyzer("StreamThingAnalyzer",
    product_to_get = cms.string('m1')
)


process.out = cms.OutputModule("EventStreamFileWriter",
    fileName = cms.untracked.string(
        'inputfile.dat'
    )
)

process.end = cms.EndPath(process.a1*process.out)

