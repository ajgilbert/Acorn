from WMCore.Configuration import Configuration
from multiprocessing import Process
config = Configuration()

PROD='vhbb-miniaod'

config.section_('General')
config.General.workArea=PROD
config.General.requestName='vhbb-miniaod-test3'

config.section_('JobType')
config.JobType.scriptExe = 'run_chain.sh'
config.JobType.psetName = 'do_nothing_cfg.py'
config.JobType.pluginName = 'PrivateMC'
config.JobType.outputFiles = ['step_3.root']
config.JobType.inputFiles = [
    'modifyCfg.py',
    'run_chain.sh',
    'vhbb-miniaod-2017/HIG-RunIIFall17DRPremix-00742_1_cfg.py',
    'vhbb-miniaod-2017/HIG-RunIIFall17DRPremix-00742_2_cfg.py',
    'vhbb-miniaod-2017/HIG-RunIIFall17MiniAODv2-00605_1_cfg.py',
    'vhbb-miniaod-2017/HIG-RunIIFall17wmLHEGS-00620_1_cfg.py',
    'chain_step_0.sh',
    'chain_step_1.sh',
    'chain_step_2.sh',
    'chain_step_3.sh'
    ]
config.JobType.disableAutomaticOutputCollection = True
config.JobType.maxMemoryMB = 6000
config.JobType.numCores = 4

config.section_('Data')
config.Data.unitsPerJob = 2000
config.Data.totalUnits = 20000
config.Data.splitting = 'EventBased'
config.Data.publication = True
#config.Data.ignoreLocality = True
config.Data.outputPrimaryDataset = PROD
config.Data.outLFNDirBase = '/store/group/cmst3/group/htautau/vhbb-testprod'
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



