#!/usr/bin/env python
import ROOT
import argparse
import re
import subprocess
import sys
import json

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

parser = argparse.ArgumentParser()

parser.add_argument('config', default='config.json',
                    help='Specifies the sample config json')
parser.add_argument('--samples', '-s', default=None,
                    help='Specifies the samples to run')
parser.add_argument('--no-cfg', action='store_true',
                    help='Do not add the preset config arguments')
parser.add_argument('--nfiles', '-n', default=1, type=int,
                    help='number of files to process')

args, passthru = parser.parse_known_args()


with open(args.config) as jsonfile:
    cfg = json.load(jsonfile)

samples = []
if args.samples is not None:
    samples = args.samples.split(',')


for s in samples:
    info = cfg['samples'][s]
    dataset = info['dataset']
    runargs = ['cmsRun'] + passthru
    if not args.no_cfg:
        runargs.extend(cfg['configs'][info['config']])
    dbs_instance = 'global'
    if 'inputDBS' in info:
        dbs_instance = info['inputDBS']
        print dbs_instance
        das_args = ['dasgoclient', '-query', 'file dataset=%s instance=prod/%s' % (dataset,  dbs_instance), '-limit', '%i' % args.nfiles]
        print ' '.join(das_args)
    infile = subprocess.check_output(das_args).strip()
    infile = ','.join([x.strip() for x in infile.split('\n')])
    runargs.append('input=%s' % infile)
    print ' '.join(runargs)
    subprocess.call(runargs)
