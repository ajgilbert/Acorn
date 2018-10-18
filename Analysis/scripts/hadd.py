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

args = parser.parse_args()
job_mgr.set_args(args)

hadddict = defaultdict(list)
for root, dirs, files in os.walk(args.dir):
    for f in files:
        basename, ext = os.path.splitext(f)
        res = re.match('^(.*)_\d+$', basename)
        if res:
            target = os.path.join(root, res.groups()[0]+ext)
            hadddict[target].append(os.path.join(root, f))
    if not args.recursive:
        break
# print hadddict

for target, inputs in hadddict.items():
    job = 'hadd -O -f6 %s %s' % (target, ' '.join(inputs))
    if args.clean:
        job += '; if [ $? == "0" ]; then rm %s; fi' % (' '.join(inputs))
    job_mgr.job_queue.append(job)

job_mgr.flush_queue()
