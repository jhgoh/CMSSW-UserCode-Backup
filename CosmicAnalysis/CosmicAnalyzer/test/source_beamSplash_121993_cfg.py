import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/data/BeamCommissioning09/Cosmics/RECO/v2/000/121/993/ECA3C304-A0D6-DE11-8BCD-0030487A18F2.root' ] );


secFiles.extend( [
       '/store/data/BeamCommissioning09/Cosmics/RAW/v1/000/121/993/D04EA868-5FD6-DE11-B372-003048D2BE08.root'] );


