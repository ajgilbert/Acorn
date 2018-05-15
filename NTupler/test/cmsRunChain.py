import argparse
import os
import json


MASTER_JOB = """#!/bin/sh
set -e
#set -x

SEED=$1

%(CONTENT)s

"""

SINGLE_JOB = """#!/bin/sh
set -e
#set -x

SEED=$1

RUN_DIR=${PWD}
echo ">> Setting RUN_DIR to ${RUN_DIR}"

cd ${CMSSW_BASE}/../
if [ ! -d chain ]; then
  mkdir chain
fi
cd chain
CHAIN_DIR=${PWD}
echo ">> Setting CHAIN_DIR to ${CHAIN_DIR}"

CMSSW_RELEASE=%(CMSSW_RELEASE)s
SCRAM_ARCH=%(SCRAM_ARCH)s

if [ "${CMSSW_RELEASE}" != "local" ]; then
    echo ">> Setting up release area for ${CMSSW_RELEASE} and arch ${SCRAM_ARCH}"
    if [ ! -d ${CMSSW_RELEASE} ]; then
      scram project CMSSW ${CMSSW_RELEASE}
    fi

    cd %(CMSSW_RELEASE)s
    eval `scramv1 runtime -sh`

    cd ${CHAIN_DIR} 
fi

%(FILE_COPY)s

%(MOD_SCRIPT)s

cmsRun -j %(JOB_REPORT)s %(RUN_SCRIPT)s

%(RETURN_FILES)s
"""

parser = argparse.ArgumentParser()
parser.add_argument('config', default='config.json',
                    help='Specify the sequence to run')
parser.add_argument('--run', action='store_true',
                    help='Actually runt the sequence')
parser.add_argument('--events', default=100,
                    help='Specify the number of events to generate')

args = parser.parse_args()

with open(args.config) as jsonfile:
    sequence = json.load(jsonfile)

master_commands = []

for i, seq in enumerate(sequence):
    release = seq['release'] 
    JOB_SCRIPT = SINGLE_JOB

    if release == 'local':
        scram_arch = os.environ['SCRAM_ARCH']
    else:
        scram_arch = seq['SCRAM_ARCH']

    run_script = seq['cfg']

    copy_arg = '\n'.join(['cp ${RUN_DIR}/%s ${CHAIN_DIR}/%s' % (x, x) for x in seq['local_files']])
    copy_arg += '\n'
    copy_arg += '\n'.join(['cp ${RUN_DIR}/%s ${CHAIN_DIR}/${CMSSW_RELEASE}/src/%s' % (x, x) for x in seq['cmssw_files']])

    events = -1
    inputFileArg = '--inputFile=step_%i.root' % (i - 1)
    outputFileArg = '--outputFile=step_%i.root' % (i)
    if i == 0:
        events = args.events
        inputFileArg = ''
    inputScript = os.path.join('${RUN_DIR}', run_script)
    outputScript = os.path.join('${CHAIN_DIR}', 'step_%i_cfg.py' % i)

    mod_arg = 'python ${RUN_DIR}/modifyCfg.py %s %s --events=%i --randomSeeds=${SEED} %s %s' % (inputScript, outputScript, events, inputFileArg, outputFileArg)

    return_files_arg = ''
    if 'output_files' in seq:
        return_files_arg = '\n'.join(['cp ${CHAIN_DIR}/%s ${RUN_DIR}/%s' % (x, x) for x in seq['output_files']])
        
    gen_script = JOB_SCRIPT % ({
        'SCRAM_ARCH': scram_arch,
        'CMSSW_RELEASE': release,
        'RUN_SCRIPT': outputScript,
        'FILE_COPY': copy_arg,
        'MOD_SCRIPT': mod_arg,
        'RETURN_FILES': return_files_arg,
        'JOB_REPORT': 'JobReport_%i.xml' % i
        })
    #print gen_script
    script_name = 'chain_step_%i.sh' % i

    with open(script_name, "w") as text_file:
        text_file.write(gen_script)

    master_commands.append('bash %s ${SEED}' % script_name)

master_script = MASTER_JOB % ({
    'CONTENT': '\n'.join(master_commands)
  })

with open('run_chain.sh', "w") as text_file:
    text_file.write(master_script)

    #os.system('chmod +x %s; ./%s' % (script_name, script_name))
