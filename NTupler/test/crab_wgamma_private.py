from WMCore.Configuration import Configuration
from multiprocessing import Process
config = Configuration()

PROD='wgamma-genonly-miniaod'

config.section_('General')
config.General.workArea=PROD
config.General.requestName='wgamma-ptA-300-inf-genonly-miniaod'

config.section_('JobType')
config.JobType.scriptExe = ''
config.JobType.psetName = ''
config.JobType.pluginName = 'PrivateMC'
config.JobType.outputFiles = []
config.JobType.inputFiles = []
config.JobType.disableAutomaticOutputCollection = True
config.JobType.maxMemoryMB = 2000
config.JobType.numCores = 1

config.section_('Data')
config.Data.unitsPerJob = 10000
config.Data.totalUnits = 5000000
config.Data.splitting = 'EventBased'
config.Data.publication = True
#config.Data.ignoreLocality = True
config.Data.outputPrimaryDataset = 'WGToMuNuG-EFT_pTA_300_inf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
config.Data.outLFNDirBase = '/store/group/cmst3/group/htautau/wgamma-testprod'
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



