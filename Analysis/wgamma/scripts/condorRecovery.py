import os
# import json
# import argparse
import subprocess
# import imp
# import re
# import pprint
import glob

indir = 'jobs/*.log'

print glob.glob(indir)

for joblog in glob.glob(indir):
    failed = []
    basename = os.path.basename(joblog)
    basename = os.path.splitext(basename)[0]
    cluster = basename[basename.rfind('.')+1:]
    basename = basename[:basename.rfind('.')]
    # print basename
    out = subprocess.check_output(['condor_history', '-userlog', joblog])
    jobresults = out.split('\n')[1:-1]
    for i, jobres in enumerate(jobresults):
        if jobres.split()[4] == 'X':
            failed.append(str(i))
            print jobres.split()
    if len(failed):
        subfile = 'jobs/condor_%s.sub' % basename
        scriptfile = 'jobs/condor_%s.sh' % basename
        rec_subfile = 'jobs-recover/condor_%s.sub' % basename
        rec_scriptfile = 'jobs-recover/condor_%s.sh' % basename
        os.system('cp %s %s' % (subfile, rec_subfile))
        os.system('cp %s %s' % (scriptfile, rec_scriptfile))
        os.system("""sed -i 's/jobs\//jobs-recover\//g' %s""" % rec_subfile)
        os.system("""sed -i 's/^queue.*$/queue arguments in (%s)/g' %s""" % ((', '.join(failed)), rec_subfile))
