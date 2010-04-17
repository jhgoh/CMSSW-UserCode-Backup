import FWCore.ParameterSet.Config as cms

dtVsRPCRecHitV = cms.EDAnalyzer("RPCPointVsRecHit",
    subDir = cms.string("RPC/RPCRecHitV/DTVsReco"),
    refHit = cms.InputTag("rpcPointProducer", "RPCDTExtrapolatedPoints"),
    recHit = cms.InputTag("rpcRecHits"),
    standAloneMode = cms.untracked.bool(True),
    rootFileName = cms.untracked.string("")
)

cscVsRPCRecHitV = cms.EDAnalyzer("RPCPointVsRecHit",
    subDir = cms.string("RPC/RPCRecHitV/CSCVsReco"),
    refHit = cms.InputTag("rpcPointProducer", "RPCCSCExtrapolatedPoints"),
    recHit = cms.InputTag("rpcRecHits"),
    standAloneMode = cms.untracked.bool(True),
    rootFileName = cms.untracked.string("")
)

trackVsRPCRecHitV = cms.EDAnalyzer("RPCPointVsRecHit",
    subDir = cms.string("RPC/RPCRecHitV/TrackVsReco"),
    refHit = cms.InputTag("rpcPointProducer", "RPCTrackExtrapolatedPoints"),
    recHit = cms.InputTag("rpcRecHits"),
    standAloneMode = cms.untracked.bool(True),
    rootFileName = cms.untracked.string("")
)

rpcPointVsRecHitValidation_step = cms.Sequence(dtVsRPCRecHitV+cscVsRPCRecHitV+trackVsRPCRecHitV)
