from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName   = 'AN-test-WWA-b'
#config.General.saveLogs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
#config.JobType.generator = 'lhe'
config.JobType.inputFiles = ['PHYS14_25_V2_All_L1FastJet_AK4PFchs.txt','PHYS14_25_V2_All_L2Relative_AK4PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK4PFchs.txt']
# Name of the CMSSW configuration file
config.JobType.psetName    = 'analysis.py'


config.section_("Data")
# This string determines the primary dataset of the newly-produced outputs.
# For instance, this dataset will be named /CrabTestSingleMu/something/USER
config.Data.inputDataset = '/WWA/qili-Q-Test-v3-7d492cb64f2cdaff326f939f96e45c96/USER'
#config.Data.inputDBS = 'global'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.totalUnits = 199
config.Data.publication = False

# This string is used to construct the output dataset name
config.Data.publishDataName = 'QiangTest-WWA-b'

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'
