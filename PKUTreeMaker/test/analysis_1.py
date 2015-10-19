import FWCore.ParameterSet.Config as cms

process = cms.Process( "TEST" )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#****************************************************************************************************#

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'PHYS14_25_V2::All'

process.load("VAJets.PKUCommon.goodMuons_cff")
process.load("VAJets.PKUCommon.goodElectrons_cff")
process.load("VAJets.PKUCommon.goodJets_cff")
process.load("VAJets.PKUCommon.goodPhotons_cff")
process.load("VAJets.PKUCommon.leptonicW_cff")

# Updates
process.goodMuons.src = "slimmedMuons"
process.goodElectrons.src = "slimmedElectrons"
process.goodAK4Jets.src = "slimmedJets"
process.goodPhotons.src = "slimmedPhotons"
process.Wtoenu.MET  = "slimmedMETs"
process.Wtomunu.MET = "slimmedMETs"

process.goodOfflinePrimaryVertex = cms.EDFilter("VertexSelector",
                                       src = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                       cut = cms.string("chi2!=0 && ndof >= 4.0 && abs(z) <= 24.0 && abs(position.Rho) <= 2.0"),
                                       filter = cms.bool(True)
                                       )

WBOSONCUT = "pt > 0.0"

process.leptonicVSelector = cms.EDFilter("CandViewSelector",
                                       src = cms.InputTag("leptonicV"),
                                       cut = cms.string( WBOSONCUT ), 
                                       filter = cms.bool(True)
                                       )

process.leptonicVFilter = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("leptonicV"),
                                       minNumber = cms.uint32(1),
                                       filter = cms.bool(True)
                                       )


process.leptonSequence = cms.Sequence(process.muSequence +
                                      process.eleSequence +
                                      process.leptonicVSequence +
                                      process.leptonicVSelector +
                                      process.leptonicVFilter )

#begin------------JEC on the fly--------
#Method 1
process.load("VAJets.PKUJets.redoPatJets_cff")
#process.goodAK4Jets.src = cms.InputTag("selectedPatJetsAK4")
process.jetSequence = cms.Sequence(process.redoPatJets+process.NJetsSequence)
#Mehod 2
jecLevelsAK4chs = [
'PHYS14_25_V2_All_L1FastJet_AK4PFchs.txt',
'PHYS14_25_V2_All_L2Relative_AK4PFchs.txt',
'PHYS14_25_V2_All_L3Absolute_AK4PFchs.txt'
]
#end------------JEC on the fly--------



#updates 2
#    process.goodOfflinePrimaryVertex.filter = False
#process.Wtomunu.cut = ''
#process.Wtoenu.cut = ''
#rocess.leptonicVSelector.filter = False
#    process.leptonicVSelector.cut = ''
#process.leptonicVFilter.minNumber = 0

print "++++++++++ CUTS ++++++++++\n"
print "Leptonic V cut = "+str(process.leptonicVSelector.cut)
print "\n++++++++++++++++++++++++++"
   
process.treeDumper = cms.EDAnalyzer("PKUTreeMaker",
                                    originalNEvents = cms.int32(1),
                                    crossSectionPb = cms.double(1),
                                    targetLumiInvPb = cms.double(1.0),
                                    PKUChannel = cms.string("VW_CHANNEL"),
                                    isGen = cms.bool(False),
                                    ak4jetsSrc = cms.string("cleanAK4Jets"),      
                                    jets = cms.string("cleanAK4Jets"),
                                    rho = cms.string("fixedGridRhoFastjetAll"),   
                                    jecAK4chsPayloadNames = cms.vstring( jecLevelsAK4chs ),
                                    photonSrc = cms.string("goodPhotons"),  
                                    leptonicVSrc = cms.string("leptonicV"),
                                    metSrc = cms.string("slimmedMETs"),
                                   # reclusteredmets = cms.InputTag("patMETs"),
                                   # pfmets = cms.InputTag("pfMet"),
                                    electronIDs = cms.InputTag("heepElectronID-HEEPV50-CSA14-25ns")
                                    electrons = cms.InputTag("slimmedElectrons"),
                                    photons = cms.InputTag("slimmedPhotons"),
                                    beamSpot = cms.InputTag("offlineBeamSpot","","RECO"),
                                    conversions = cms.InputTag("reducedEgamma","reducedConversions","PAT"),

                                    )


process.analysis = cms.Path(
                            process.goodOfflinePrimaryVertex +
                            process.leptonSequence +
                            process.photonSequence +
                            process.jetSequence +
                            process.treeDumper)

### Source
process.load("VAJets.PKUCommon.data.RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8")
process.source.fileNames = ["/store/user/qili/WWA/Q-Test-v3/150318_075359/0000/JME-Phys14DR-00001_MINIAOD_99.root","/store/user/qili/WWA/Q-Test-v3/150318_075359/0000/JME-Phys14DR-00001_MINIAOD_98.root","/store/user/qili/WWA/Q-Test-v3/150318_075359/0000/JME-Phys14DR-00001_MINIAOD_97.root","/store/user/qili/WWA/Q-Test-v3/150318_075359/0000/JME-Phys14DR-00001_MINIAOD_96.root"]

process.maxEvents.input = -1

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.MessageLogger.cerr.FwkReport.limit = 99999999

process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("treePKU_TT_xwh_7.root")
                                    fileName = cms.string("treePKU.root")
#                                   fileName = cms.string("treePKU_MWp_3000_bb_xwh_3_noGen_cuts_1.root")
#                                   fileName = cms.string("treePKU_MWp_4000_bb_xwh.root")
#                                   fileName = cms.string("WJetsToLNu_13TeV-madgraph-pythia8-tauola_all.root")
#                                   fileName = cms.string("JME-Fall13-00001_py8_AN_type2_allcuts_DIYJOB.root")
#                                   fileName = cms.string("treePKU.root")
                                   )
