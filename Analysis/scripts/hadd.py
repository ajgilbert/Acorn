from Acorn.Analysis.jobs import Jobs
import argparse
import os
import re
from collections import defaultdict
# import json

job_mgr = Jobs()
parser = argparse.ArgumentParser()
job_mgr.attach_job_args(parser)

parser.add_argument('dir',
                    help='root directory')
parser.add_argument('--recursive', '-r', action='store_true',
                    help='hadd recursively in subdirectories')
parser.add_argument('--clean', '-c', action='store_true',
                    help='Delete input files if hadd was successful')
parser.add_argument('--remote-dir', default=None,
                    help='Directory prefix for remote read/copy')

# root://eoscms.cern.ch//store/cmst3/user/agilbert/190205/wgamma_2016_v3
# /eos/cms/store/cmst3/user/agilbert/190205 wgamma_2016_v3/WGamma_WGToLNuG-amcatnloFXFX.root
args = parser.parse_args()
job_mgr.set_args(args)

indir = args.dir
outdir = args.dir
is_remote = False
if args.remote_dir is not None:
    is_remote = True
    outdir = args.remote_dir

hadddict = defaultdict(list)
cleandict = defaultdict(list)
for root, dirs, files in os.walk(indir):
    for f in files:
        basename, ext = os.path.splitext(f)
        res = re.match('^(.*)_\d+$', basename)
        if res:
            target_dir = root
            if is_remote:
                target_dir = target_dir.replace(indir, outdir)
            target = os.path.join(target_dir, res.groups()[0]+ext)
            hadddict[target].append(os.path.join(target_dir, f))
            cleandict[target].append(os.path.join(root, f))
    if not args.recursive:
        break
# print hadddict

for target, inputs in hadddict.items():
    actual_target = target
    if is_remote:
        actual_target = '$TMPDIR/%s' % os.path.basename(target)
    job = 'hadd -O -f6 %s %s' % (actual_target, ' '.join(inputs))
    if is_remote:
        job += ' && xrdcp %s %s && rm %s' % (actual_target, target, actual_target)
    if args.clean:
        job += ' && rm %s' % (' '.join(cleandict[target]))
    job_mgr.job_queue.append(job)

job_mgr.flush_queue()
