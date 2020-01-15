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
parser.add_argument('--tmp', '-t', default='',
                    help='Write the job output in a temporary location')
parser.add_argument('--merge-tasks', '-m', action='store_true',
                    help='Merge all tasks into one big submission')

args = parser.parse_args()
job_mgr.set_args(args)

with open(args.config) as jsonfile:
    incfg = json.load(jsonfile)

full_outdir = os.path.join(args.outdir, args.production)
sequences = args.sequences.split(',')

if not full_outdir.startswith('root://'):
    print '>> Creating directory %s' % full_outdir
    os.system('mkdir -p %s' % full_outdir)

SAMPLES = incfg['samples']

task = job_mgr.task_name

for sample in SAMPLES:
    if args.attributes not in SAMPLES[sample]['attributes']:
        continue
    filtered_seqs = list(sequences)
    if 'sequences' in SAMPLES[sample]:
        filtered_seqs = [X for X in filtered_seqs if X in SAMPLES[sample]['sequences']]
    if len(filtered_seqs) == 0:
        continue
    cfg = {
        'sequences': filtered_seqs,
        'filelists': ['filelists/%s_%s.txt' % (args.production, x) for x in SAMPLES[sample]['inputs']],
        'outdir': full_outdir,
        'attributes': SAMPLES[sample]['attributes']
    }
    output_cfgs = []
    for seq in filtered_seqs:
        cfg['output_%s' % seq] = '%s_%s.root' % (seq, sample)
        output_cfgs.append('output_%s' % seq)
    cfg.update(incfg['config'])
    if 'stitching' in SAMPLES[sample]:
        cfg['stitching'] = SAMPLES[sample]['stitching']

    files_per_job = args.files_per_job
    if 'data' in SAMPLES[sample]['attributes']:
        files_per_job *= 2
    job_mgr.add_filelist_split_jobs(
        prog='Acorn_WGAnalysis',
        cfg=cfg,
        files_per_job=files_per_job,
        output_cfgs=output_cfgs,
        tempdir=args.tmp)

    if not args.merge_tasks:
        job_mgr.task_name = task + '-' + sample
        job_mgr.flush_queue()

if args.merge_tasks:
    job_mgr.flush_queue()