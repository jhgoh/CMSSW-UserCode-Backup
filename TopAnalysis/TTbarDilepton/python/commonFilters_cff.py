import FWCore.ParameterSet.Config as cms

from CommonTools.RecoAlgos.HBHENoiseFilter_cfi import *

noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25),
)

goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter", 
    src = cms.InputTag('offlinePrimaryVertices'),
    filterParams =  cms.PSet(
        minNdof = cms.double(4.),
        maxZ    = cms.double(24.), 
        maxRho  = cms.double(2.)
    ),
    filter = cms.bool(True),
)

commonSequenceForData = cms.Sequence(
    goodOfflinePrimaryVertices
#  + noscraping ## Noscrapping does not work in cmgTuple since no generalTracks in the event
)

commonSequenceForMC = cms.Sequence(
    goodOfflinePrimaryVertices
)
