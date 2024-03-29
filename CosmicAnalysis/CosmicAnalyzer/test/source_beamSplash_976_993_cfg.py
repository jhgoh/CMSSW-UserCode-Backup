import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121981.0001.Error.storageManager.03.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121981.0001.Error.storageManager.01.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121981.0001.Error.storageManager.02.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121981.0001.Error.storageManager.04.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.02.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.00.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.02.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121982.0001.Error.storageManager.01.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.06.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.00.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.03.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.01.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.07.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.05.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.05.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.05.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.06.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.04.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.07.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.03.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.02.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.01.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.06.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.04.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.02.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.03.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.00.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.06.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.04.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.04.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.02.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.02.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.00.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.05.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.06.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.00.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.01.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.04.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.07.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.01.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.03.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.07.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.05.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.01.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121976.0001.Error.storageManager.07.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.03.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/ALCA_recoverreco_Data.00121980.0001.Error.storageManager.01.0002.dat.root',


'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/data/Express/121/943/DQM_V0001_R000121943__StreamExpress__BeamCommissioning09-Express-v2__DQM.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/data/Express/121/949/DQM_V0001_R000121949__StreamExpress__BeamCommissioning09-Express-v2__DQM.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/data/Express/121/960/DQM_V0001_R000121960__StreamExpress__BeamCommissioning09-Express-v2__DQM.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/data/Express/121/962/DQM_V0001_R000121962__StreamExpress__BeamCommissioning09-Express-v2__DQM.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/data/Express/121/964/DQM_V0001_R000121964__StreamExpress__BeamCommissioning09-Express-v2__DQM.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/data/Express/121/993/DQM_V0001_R000121993__StreamExpress__BeamCommissioning09-Express-v2__DQM.root'

#       '/store/express/BeamCommissioning09/StreamExpress/ALCARECO/v2/000/121/943/0A8735F7-08D6-DE11-8240-000423D99B3E.root',
#       '/store/express/BeamCommissioning09/StreamExpress/ALCARECO/v2/000/121/943/047A1CC1-0BD6-DE11-A910-001617C3B76A.root'
	 ] );


secFiles.extend( [
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121981.0001.Error.storageManager.03.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121981.0001.Error.storageManager.01.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121981.0001.Error.storageManager.02.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121981.0001.Error.storageManager.04.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.02.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.00.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.02.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121982.0001.Error.storageManager.01.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.06.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.00.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.03.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.01.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.07.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.05.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.05.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.05.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.06.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.04.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.07.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.03.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.02.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.01.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.06.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.04.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.02.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.03.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.00.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.06.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.04.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.04.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.02.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.02.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.00.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.05.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.06.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.00.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.01.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.04.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.07.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.01.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.03.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.01.0002.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.05.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.07.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121976.0001.Error.storageManager.07.0001.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.03.0000.dat.root',
'file:/afs/cern.ch/cms/CAF/CMSCOMM/COMM_GLOBAL/ErrorRecover/recoverreco_Data.00121980.0001.Error.storageManager.01.0002.dat.root'


	 ] );


