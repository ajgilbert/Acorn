from Acorn.Analysis.jobs import Jobs
import argparse
import os
import json

job_mgr = Jobs()
parser = argparse.ArgumentParser()
job_mgr.attach_job_args(parser)

parser.add_argument('config',
                    help='Config file')
parser.add_argument('--files-per-job', '-f', type=int, default=100,
                    help='Number of files to process per job')
parser.add_argument('--attributes', '-a', default='default',
                    help='Only run the samples matching the chosen attribute')
parser.add_argument('--outdir', '-o', default='output',
                    help='Main output directory')
parser.add_argument('--production', '-p', default='wgamma_2018_v1',
                    help='Production tag, used as filelist prefix and output path')
parser.add_argument('--sequences', '-s', default='WGamma',
                    help='List of sequences to run')

args = parser.parse_args()
job_mgr.set_args(args)

with open(args.config) as jsonfile:
    incfg = json.load(jsonfile)

full_outdir = os.path.join(args.outdir, args.production)
sequences = args.sequences.split(',')

for seq in sequences:
    makedir = os.path.join(full_outdir, seq)
    print '>> Creating directory %s' % makedir
    os.system('mkdir -p %s' % makedir)

SAMPLES = incfg['samples']

task = job_mgr.task_name

for sample in SAMPLES:
    if args.attributes not in SAMPLES[sample]['attributes']:
        continue
    cfg = {
        'sequences': sequences,
        'filelists': ['filelists/%s_%s.txt' % (args.production, x) for x in SAMPLES[sample]['inputs']],
        'output': '%s.root' % sample,
        'outdir': full_outdir,
        'attributes': SAMPLES[sample]['attributes']
    }
    cfg.update(incfg['config'])

    job_mgr.add_filelist_split_jobs(
        prog='Acorn_WGAnalysis',
        cfg=cfg,
        files_per_job=args.files_per_job,
        output_cfgs=['output'])
    job_mgr.task_name = task + '-' + sample
    job_mgr.flush_queue()
