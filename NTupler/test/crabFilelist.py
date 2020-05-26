import argparse
import json
import sys
# from WMCore.Configuration import Configuration
from CRABAPI.RawCommand import crabCommand
from httplib import HTTPException
import subprocess
import os
import re

parser = argparse.ArgumentParser()


parser.add_argument('dir',
                    help='Specifies the crab directory to process')
parser.add_argument('--replace-xrootd', default='eoscms.cern.ch',
                    help='Replace the global redirector with the actual xrootd server')
parser.add_argument('--prefix', default='',
                    help='Add this prefix to the output filelist')
parser.add_argument('--filename', default='EventTree.root',
                    help='Output filename')

args = parser.parse_args()

taskname = args.dir.replace('crab_','')
filelist_name = '%s%s.txt' % (args.prefix, taskname)
if os.path.isfile(filelist_name):
    print '>> Filelist %s already exists' % filelist_name
    sys.exit(0)
print '>> Processing directory %s' % args.dir
res = crabCommand('status', dir=args.dir)
done_ids = set()
all_ids = set()
file_ids = set()
for status, index in res['jobList']:
    all_ids.add(index)
    if status == 'finished':
        done_ids.add(index)
print '>> Jobs finished: %i/%i' % (len(done_ids), len(all_ids))

if len(done_ids) == 0:
    print '>> No jobs have finished yet, quitting'
    sys.exit(1)

test_job = list(done_ids)[0]

print '>> Getting xrootd path for a random job...'
res = crabCommand('getoutput', dir=args.dir, dump=True, xrootd=True, jobids=test_job)
xrootd_dir = '/'.join(res['xrootd'][0].split('/')[:-2])
xrootd_dir = xrootd_dir.replace('cms-xrd-global.cern.ch', args.replace_xrootd)

print '>> The main xrootd directory is: %s' % xrootd_dir

filename = args.filename 
filename_split = os.path.splitext(filename) 
filename_re = re.compile(r'%s_(\d+)\%s' % filename_split)

print '>> Getting full PFN path for a random job...'
res = crabCommand('getoutput', dir=args.dir, dump=True, jobids=test_job)
pfn = res['pfn']
main_dir = '/'.join(pfn[0].split('/')[:-2])

print '>> The main PFN directory is: %s' % main_dir

print '>> Looking for %s files in this directory...' % args.filename
# Workaround for running under CC7 - the cvmfs/CMSSW python version breaks gfal commands
# We have to hide the contents of LD_LIBRARY_PATH:
env_for_gfal = os.environ.copy()
env_for_gfal['LD_LIBRARY_PATH'] = ""

out = subprocess.check_output(['gfal-ls', main_dir], env=env_for_gfal)
subdirs = [x for x in out.split('\n') if x != '']
all_matched_files = []
for subdir in subdirs:
    out = subprocess.check_output(['gfal-ls', '/'.join([main_dir, subdir])], env=env_for_gfal)
    files = [x for x in out.split('\n') if x != '']
    for f in files:
        res = re.match(filename_re, f)
        if re.match(filename_re, f):
            index = res.groups()[0]
            if index not in done_ids:
                print '>> Found output file for job index %i, not yet in finished state - skipping\n' % int(index)
                continue
            file_ids.add(index)
            all_matched_files.append('/'.join([xrootd_dir, subdir, f]))

missing_files = done_ids - file_ids

if len(missing_files):
    print '>> ERROR: the following jobids are reported as finished by crab, but no output file was found'
    sys.exit(1)
else:
    print '>> All %i files were found, writing filelist %s' % (len(done_ids), filelist_name) 

with open(filelist_name, "w") as text_file:
    text_file.write('\n'.join(all_matched_files))


