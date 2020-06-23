import argparse
import json
import sys
from CRABClient.UserUtilities import config as Configuration
from multiprocessing import Process
from CRABAPI.RawCommand import crabCommand
from httplib import HTTPException

parser = argparse.ArgumentParser()

parser.add_argument('config', default='config.json',
                    help='Specifies the samples to run')
parser.add_argument('--submit', action='store_true',
                    help='Actually do the submission')
parser.add_argument('--unitsPerJob', '-u', default=200000, type=int,
                    help='Specify the number of units per job')
parser.add_argument('--psetName', '-p', default='test/wgamma_cfg.py',
                    help='PSet to use')
parser.add_argument('--label', '-l', default='production',
                    help='Production label')
parser.add_argument('--verbosity', '-v', type=int, default=0,
                    help='Level of verbosity')
parser.add_argument('--attribute', '-a', default='submit',
                    help='Submit samples having this attribute')
args = parser.parse_args()

config = Configuration()

with open(args.config) as jsonfile:
    samples = json.load(jsonfile)

crab_settings = samples['crabSettings']

config.section_('General')
config.General.workArea = args.label
#config.General.requestName = ''


config.section_('JobType')
config.JobType.psetName = args.psetName
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['EventTree.root']
config.JobType.pyCfgParams = []
config.JobType.allowUndistributedCMSSW = True
if 'maxMemoryMB' in crab_settings:
    config.JobType.maxMemoryMB = int(crab_settings['maxMemoryMB'])
if 'numCores' in crab_settings:
    config.JobType.numCores = int(crab_settings['numCores'])
if 'maxJobRuntimeMin' in crab_settings:
    config.JobType.maxJobRuntimeMin = int(crab_settings['maxJobRuntimeMin'])
if 'inputFiles' in crab_settings:
    config.JobType.inputFiles = [str(x) for x in crab_settings['inputFiles']]

config.section_('Data')
#config.Data.inputDataset = ''
config.Data.unitsPerJob = args.unitsPerJob
config.Data.splitting = 'EventAwareLumiBased'
config.Data.publication = False
config.Data.ignoreLocality = False
if 'outputDir' in crab_settings:
  config.Data.outLFNDirBase = '%s%s' % (crab_settings['outputDir'],args.label)
if 'allowNonValidInputDataset' in crab_settings:
  config.Data.allowNonValidInputDataset = bool(crab_settings['allowNonValidInputDataset'])

config.section_('User')
if 'voGroup' in crab_settings :
  config.User.voGroup = crab_settings['voGroup']

config.section_('Site')
config.Site.storageSite = crab_settings['storageSite']
if 'whiteList' in crab_settings:
  config.Site.whitelist = [str(x) for x in crab_settings['whiteList']]
if 'blacklist' in crab_settings:
  config.Site.blacklist = [str(x) for x in crab_settings['blacklist']]
#config.Site.whitelist = ['T2_UK_London_IC', 'T2_CH_CERN', 'T2_FR_GRIF_LLR', 'T2_UK_SGrid_Bristol', 'T3_US_FNALLPC', 'T2_DE_DESY', 'T2_IT_Bari', 'T2_BE_IIHE', 'T2_US_UCSD', 'T2_US_MIT', 'T2_IT_Pisa', 'T2_US_Wisconsin', 'T2_US_Florida', 'T2_IT_Rome','T2_FR_IPHC']


def submit(config):
    try:
        crabCommand('submit', config=config)
    except HTTPException, hte:
        print hte.headers
        sys.exit(1)

print '>> Production with label: %s' % args.label

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################
for sample in sorted(samples['samples']):
    info = samples['samples'][sample]
    if args.attribute in info['attributes']:
        print '>> Will submit: %-40s: %s' % (sample, info['dataset'])
    else:
        if args.verbosity >= 1:
            print '>> Skipping: %-40s: %s' % (sample, info['dataset'])
        continue
    # Want to make a copy of the configuration to modify freely
    # Apparently this is the way to do it...
    config2 = Configuration()
    config2 += config

    config2.General.requestName = '%s_%s' % (args.label, sample)
    config2.Data.outputDatasetTag = sample
    config2.Data.inputDataset = info['dataset']
    config2.JobType.pyCfgParams = [str(x) for x in samples['configs'][info['config']]]
    if 'runRange' in info:
        config2.Data.runRange = str(info['runRange'])
    if 'lumiMask' in info:
        config2.Data.lumiMask = str(info['lumiMask'])
    if 'useFileSplitting' in info:
        config2.Data.splitting = 'FileBased'
        config2.Data.unitsPerJob = info['useFileSplitting']
        # config2.General.instance = 'preprod' # temporary - see https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/3623/1/1.html
    if 'inputDBS' in info:
        config2.Data.inputDBS = str(info['inputDBS'])
    if args.verbosity >= 1:
        print config2
    if (args.submit):
        p = Process(target=submit, args=(config2,))
        p.start()
        p.join()
        # submit(config2)
