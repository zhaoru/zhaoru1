# or get from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
import FWCore.ParameterSet.Config as cms
process = cms.Process("jectxt")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration.StandardSequences.Services_cff')
# define your favorite global tag
process.GlobalTag.globaltag = 'PHYS14_25_V2::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.source = cms.Source("EmptySource")
process.readAK4PF = cms.EDAnalyzer('JetCorrectorDBReader',
# below is the communication to the database
payloadName = cms.untracked.string('AK4PFchs'),
# this is used ONLY for the name of the printed txt files. You can use any name that you like,
# but it is recommended to use the GT name that you retrieved the files from.
globalTag = cms.untracked.string('PHYS14_25_V2::All'),
printScreen = cms.untracked.bool(False),
createTextFile = cms.untracked.bool(True)
)
process.p = cms.Path(process.readAK4PF )
