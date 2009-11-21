import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
	        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/121/943/3C83AB0D-12D6-DE11-9280-001D09F2525D.root'
	] );


secFiles.extend( [
       '/store/data/BeamCommissioning09/MinimumBias/RAW/v1/000/121/943/0841DB21-0DD6-DE11-8B4F-001D09F232B9.root'

	 ] );


