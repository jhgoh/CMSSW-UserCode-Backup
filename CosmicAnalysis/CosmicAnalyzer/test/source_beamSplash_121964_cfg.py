import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/data/BeamCommissioning09/Cosmics/RECO/v2/000/121/964/10056149-43D6-DE11-BDAF-0030487C5CFA.root' ] );


secFiles.extend( [
       '/store/data/BeamCommissioning09/Cosmics/RAW/v1/000/121/964/5610AF0C-31D6-DE11-A9F8-001617C3B6DE.root'] );


