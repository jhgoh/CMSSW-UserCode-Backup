import FWCore.ParameterSet.Config as cms

source = cms.Source("EmptySource")

hppMass = '1000.0'

from Configuration.Generator.PythiaUESettings_cfi import *
generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.),
    comEnergy = cms.double(10000.0),
#    crossSection = cms.untracked.double(7.1),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        processParameters = cms.vstring(
            'PMAS(353,1)=%s          !mass of H_L' % hppMass, 
            'PMAS(354,1)=%s          !mass of H_R' % hppMass, 
            'MSEL=0                  !(D=1) to select between full user control (0, then use MSUB) and some preprogrammed alternative: QCD hight pT processes (1, then ISUB=11, 12, 13, 28, 53, 68), QCD low pT processes (2, then ISUB=11, 12, 13, 28, 53, 68, 91, 92, 94, 95)', 
            'MSUB(349)=1             !H_L H_L', 
            'MSUB(350)=1             !H_R H_R', 
            'MSUB(351)=0             !H_L ff (WW fusion)', 
            'MSUB(352)=0             !H_R ff (WW fusion)', 
            'CKIN(45)=5.             !high mass cut on m2 in 2 to 2 process Registered by Chris.Seez@cern.ch', 
            'MSTP(25)=2              !Angular decay correlations in H->ZZ->4fermions Registered by Alexandre.Nikitenko@cern.ch', 
            'PARP(181)= 0.0          !H++ coupling with ee', 
            'PARP(182)= 0.1          !H++ coupling with emu', 
            'PARP(183)= 0.0          !H++ coupling with etau', 
            'PARP(184)= 0.0          !H++ coupling with mue', 
            'PARP(185)= 0.0          !H++ coupling with mumu', 
            'PARP(186)= 0.0          !H++ coupling with mutau', 
            'PARP(187)= 0.0          !H++ coupling with taue', 
            'PARP(188)= 0.0          !H++ coupling with taumu', 
            'PARP(189)= 0.0          !H++ coupling with taumutau', 
            'PARP(190)= 0.0          !H++_L coupling with W', 
            'PARP(191)= 0.0          !H++_R coupling with W'
         ),
        # This is a vector of ParameterSet names to be read, in this order
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)

ProductionFilterSequence = cms.Sequence(generator)

configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.2 $'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/UserCode/JHGoh/HiggsAnalysis/DoublyChargedHiggs/python/PYTHIA6_HppMuMu_M1000_10TeV_cff.py,v $'),
    annotation = cms.untracked.string('PYTHIA6-H++ to MuMu at 10TeV')
)


