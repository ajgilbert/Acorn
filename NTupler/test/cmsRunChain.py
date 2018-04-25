import argparse
import os

EXTERNAL_JOB_PREFIX = """#!/bin/sh
set -e
set -x
cd %(CHAIN_BASE)s
mkdir -p chain
cd chain
export SCRAM_ARCH=%(SCRAM_ARCH)s

if [ ! -d %(CMSSW_RELEASE)s ]; then
  scram project CMSSW %(CMSSW_RELEASE)s
fi

cd %(CMSSW_RELEASE)s
eval `scramv1 runtime -sh`

cd %(RUN_BASE)s
cmsRun %(RUN_SCRIPT)s %(SCRIPT_ARGS)s
"""

LOCAL_JOB_PREFIX = """#!/bin/sh
set -e
set -x
cmsRun %(RUN_SCRIPT)s %(SCRIPT_ARGS)s
"""



run_base = os.environ['PWD']
cmssw_base = os.environ['CMSSW_BASE']
chain_base = os.path.dirname(cmssw_base)

parser = argparse.ArgumentParser()
parser.add_argument('sequence', nargs='*',
                    help='Specify the sequence to run')

args = parser.parse_args()

for i, seq in enumerate(args.sequence):
    all_settings = seq.split(':')
    release_settings = all_settings[0].split(',')
    cmssw_release = release_settings[0]
    JOB_PREFIX = EXTERNAL_JOB_PREFIX

    if cmssw_release == 'local':
        JOB_PREFIX = LOCAL_JOB_PREFIX
    else:
        scram_arch = release_settings[1]

    run_script = all_settings[1]
    if len(all_settings) > 2:
        script_args = ' '.join(all_settings[2:])
    else:
        script_args = ''

    gen_script = JOB_PREFIX % ({
        'CHAIN_BASE': chain_base,
        'SCRAM_ARCH': scram_arch,
        'CMSSW_RELEASE': cmssw_release,
        'RUN_BASE': run_base,
        'RUN_SCRIPT': run_script,
        'SCRIPT_ARGS': script_args
        })
    print gen_script
    script_name = 'chain_step_%i.sh' % i

    with open(script_name, "w") as text_file:
        text_file.write(gen_script)

    os.system('chmod +x %s; ./%s' % (script_name, script_name))

