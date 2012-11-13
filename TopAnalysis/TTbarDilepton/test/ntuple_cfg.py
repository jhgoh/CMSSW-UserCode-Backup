import FWCore.ParameterSet.Config as cms
import sys, os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
#process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.GlobalTag.globaltag = "START52_V12::All"

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
from TopAnalysis.TTbarDilepton.dataset_cff import *
dataset, section, nFiles = parseJobSectionOption()
files = loadDataset(dataset)
begin, end = calculateRange(files, section, nFiles)
process.source.fileNames = files[begin:end]

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple/unmerged/ntuple_%s_%03d.root" % (dataset, section)),
)

process.genParticleCount = cms.EDFilter("GenParticleCountFilter",
    src = cms.InputTag("genParticlesPruned"),
    cut = cms.string("status == 3"),
    ids = cms.vint32(13, -13, 11, -11),
    minNumber = cms.uint32(1), 
    maxNumber = cms.uint32(1),
)

process.genParticleTauVeto =cms.EDFilter("GenParticleCountFilter",
    src = cms.InputTag("genParticlesPruned"),
    cut = cms.string("status == 3"),
    ids = cms.vint32(15, -15),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(0),
)

process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
    src = cms.InputTag("genParticlesPruned"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False)
)

process.load("Configuration.StandardSequences.Generator_cff")
process.genParticlesForJetsNoNu.src = "genParticlesPruned"

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltHighLevel.throw = False
process.hltEE = process.hltHighLevel.clone()
process.hltMM = process.hltHighLevel.clone()
process.hltME = process.hltHighLevel.clone()

process.load("TopAnalysis.GeneratorTools.genJetAssociation_cff")
process.load("TopAnalysis.GeneratorTools.lumiWeight_cff")
process.load("TopAnalysis.TTbarDilepton.commonFilters_cff")
process.load("TopAnalysis.TTbarDilepton.eventTupleProducer_cfi")

process.mm = process.eventTuple.clone()
process.ee = process.eventTuple.clone()
process.me = process.eventTuple.clone()

process.mm.muon.minNumber = 2
process.ee.electron.minNumber = 2
process.me.muon.minNumber = 1
process.me.electron.minNumber = 1

production, primaryDS = dataset.split('-')
from TopAnalysis.TTbarDilepton.trigger_cff import *

if isRealData(dataset):
    process.commonSequence = cms.Sequence(
        process.commonSequenceForData
    )

    if 'DoubleMu' == primaryDS:
        process.hltMM.HLTPaths = HLTPaths[dataset]
        process.p = cms.Path(process.commonSequence + process.hltMM + process.mm)
    elif 'DoubleElectron' == primaryDS:
        process.hltEE.HLTPaths = HLTPaths[dataset]
        process.p = cms.Path(process.commonSequence + process.hltEE + process.ee)
    elif 'MuEG' == primaryDS:
        process.hltME.HLTPaths = HLTPaths[dataset]
        process.p = cms.Path(process.commonSequence + process.hltME + process.me)

else:
    process.mm.doMCMatch = True
    process.ee.doMCMatch = True
    process.me.doMCMatch = True

    process.hltMM.HLTPaths = HLTPaths["%s-DoubleMu"       % production]
    process.hltEE.HLTPaths = HLTPaths["%s-DoubleElectron" % production]
    process.hltME.HLTPaths = HLTPaths["%s-MuEG"           % production]

    process.commonSequence = cms.Sequence(
        process.commonSequenceForMC
    #  + process.genParticleCount + process.genParticleTauVeto
      + process.genParticlesForJetsNoNu * process.ak5GenJetsNoNu
      + process.recoToGenJet + process.genJetToPartons
      + process.lumiWeight
    )

    prod, sample = dataset.split('-')

    process.pMM = cms.Path(process.commonSequence + process.hltMM + process.mm)
    process.pEE = cms.Path(process.commonSequence + process.hltEE + process.ee)
    process.pME = cms.Path(process.commonSequence + process.hltME + process.me)

