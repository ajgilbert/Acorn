from WMCore.Configuration import Configuration
from multiprocessing import Process
config = Configuration()

PROD='wgamma-testprod-miniaod-2017'

config.section_('General')
config.General.workArea=PROD
config.General.requestName='wgamma-miniaod-2017'

config.section_('JobType')
config.JobType.scriptExe = ''
config.JobType.psetName = ''
config.JobType.pluginName = 'PrivateMC'
config.JobType.outputFiles = []
config.JobType.inputFiles = []
config.JobType.disableAutomaticOutputCollection = True
config.JobType.maxMemoryMB = 6000
config.JobType.numCores = 4

config.section_('Data')
config.Data.unitsPerJob = 5000
config.Data.totalUnits = 5000000
config.Data.splitting = 'EventBased'
config.Data.publication = True
#config.Data.ignoreLocality = True
config.Data.outputPrimaryDataset = 'WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8'
config.Data.outLFNDirBase = '/store/group/cmst3/group/htautau/wgamma-miniaod-2017'
#config.Data.inputDBS = 'phys03'

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_CH_CERN', 'T1_US_FNAL', 'T2_DE_DESY', 'T1_DE_KIT']
#config.Site.blacklist = ['T2_US_Purdue']

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).

    def submit(config):
        #print config
        #return
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print hte.headers

    p = Process(target=submit, args=(config,))
    p.start()
    p.join()



