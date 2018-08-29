import os
import json
import argparse
import subprocess
import imp
import re
import pprint
import FWCore.ParameterSet.Config as cms

def GetConfigForRun(run):
    # hltConfigFromDB --runNumber 276811 --nomodules --nopaths --noes --noservices --nooutput --nopsets
    out = subprocess.check_output(['hltConfigFromDB', '--runNumber', str(run), '--nomodules', '--nopaths', '--noes', '--nooutput', '--nopsets'])
    # for some reason the last line isn't valid config, so I'll remove it
    run_config = {}
    exec(out, run_config)
    return run_config['process']

def GetListOfPrescaledPaths(process):
    result = []
    lvl1Labels = [x for x in process.PrescaleService.lvl1Labels]
    for pset in process.PrescaleService.prescaleTable:
        prescales = [x for x in pset.prescales]
        unprescaled = True
        for label, ps in zip(lvl1Labels, prescales):
            if label not in ['Emergency'] and ps != 1:
                unprescaled = False
                break
        if not unprescaled:
            result.append(pset.pathName)
    return result


parser = argparse.ArgumentParser()
# parser.add_argument('input')
parser.add_argument('--lumi-json', '-l', default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-316271_13TeV_PromptReco_Collisions18_JSON.txt')
parser.add_argument('--filter-streams', default='^Physics.*')
# parser.add_argument('--event', '-e', default=1, type=int)
# parser.add_argument('--region', '-r', default='p')
# parser.add_argument('--window', default=0.1, type=float)

# parser.add_argument('--bkg-model', default='Exponential')
# parser.add_argument('--title', default='Muon ID Efficiency')
# parser.add_argument('--postfix', default='')
# parser.add_argument('--plot-dir', '-p', default='./')
# parser.add_argument('--bin-replace', default=None) #(100,2.3,80,2.3)
args = parser.parse_args()


lumi_json = {}
if args.lumi_json is not None:
    with open(args.lumi_json) as jsonfile:
      lumi_json = json.load(jsonfile)

lumi_list = sorted([int(x) for x in lumi_json.keys()])
print lumi_list

for i, run in enumerate(lumi_list):
    print '>> Analysing HLT config for run %i' % run
    if i == 0 or i == (len(lumi_list) - 1):
        process = GetConfigForRun(run)
        prescaled_paths = GetListOfPrescaledPaths(process)
        for stream in process.streams.parameterNames_():
          if re.match(args.filter_streams, stream):
            print '>> Found matching stream: %s' % stream
            print '>> Containing datasets:'
            datasets = [x for x in getattr(process.streams, stream)]
            pprint.pprint(datasets, indent=2)
            for dataset in datasets:
                paths = [x for x in getattr(process.datasets, dataset)]
                prescaled = [bool(x in prescaled_paths) for x in paths]
                print '>> Dataset %s contains [paths, prescaled]:' % dataset
                for path, ps in zip(paths, prescaled):
                    print '  %-80s %i' % (path, ps)


"""
paths_2016_MC = [
    'HLT_IsoMu22_v',
    'HLT_IsoMu22_eta2p1_v',
    'HLT_IsoTkMu22_v',
    'HLT_IsoTkMu22_eta2p1_v',
    ]

for path in paths_2016_MC:
    print 'brilcalc trg --prescale --hltpath "%s*" -o %s.csv' % (path, path)
    os.system('brilcalc trg --prescale --hltpath "%s*" -o %s.csv' % (path, path))
    """
