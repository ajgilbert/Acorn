import argparse
import json
import sys
from WMCore.Configuration import Configuration
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


config.section_('General')
config.General.workArea = args.label
#config.General.requestName = ''


config.section_('JobType')
config.JobType.psetName = args.psetName
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['EventTree.root']
config.JobType.pyCfgParams = []


config.section_('Data')
#config.Data.inputDataset = ''
config.Data.unitsPerJob = args.unitsPerJob
config.Data.splitting = 'EventAwareLumiBased'
config.Data.publication = False
config.Data.ignoreLocality = False
config.Data.outLFNDirBase = '/store/group/cmst3/group/htautau/%s' % args.label


config.section_('User')
# config.User.voGroup = 'dcms'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_UK_London_IC', 'T2_CH_CERN', 'T2_FR_GRIF_LLR', 'T2_UK_SGrid_Bristol', 'T3_US_FNALLPC', 'T2_DE_DESY', 'T2_IT_Bari', 'T2_BE_IIHE', 'T2_US_UCSD', 'T2_US_MIT', 'T2_IT_Pisa', 'T2_US_Wisconsin', 'T2_US_Florida', 'T2_IT_Rome','T2_FR_IPHC']


def submit(config):
    try:
        crabCommand('submit', config=config)
    except HTTPException, hte:
        print hte.headers
        sys.exit(1)

print '>> Production with label: %s' % args.label

with open(args.config) as jsonfile:
    samples = json.load(jsonfile)
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
    config.General.requestName = sample
    config.Data.inputDataset = info['dataset']
    config.JobType.pyCfgParams = [str(x) for x in samples['configs'][info['config']]]
    if 'runRange' in info:
        config.Data.runRange = str(info['runRange'])
    if args.verbosity >= 1:
        print config
    if (args.submit):
        submit(config)