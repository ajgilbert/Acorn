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

DO_NOTHING_CFG="""import FWCore.ParameterSet.Config as cms
process = cms.Process("MAIN")

process.source = cms.Source("EmptySource")
process.options = cms.untracked.PSet(

)

process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('file:%(FILENAME)s')
)
%(OTHER_SETTINGS)s

process.output_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.output_step)
"""

parser = argparse.ArgumentParser()
parser.add_argument('config', default='config.json',
                    help='Specify the sequence to run')
parser.add_argument('--run', action='store_true',
                    help='Actually runt the sequence')
parser.add_argument('--events', default=100, type=int,
                    help='Specify the number of events to generate')
parser.add_argument('--crab', action='store_true',
                    help='Prepare for running with crab')
parser.add_argument('--crab-cfg', default=None,
                    help='Template crab config')
parser.add_argument('--label', default='test',
                    help='Prepare for running with crab')

args = parser.parse_args()

with open(args.config) as jsonfile:
    sequence = json.load(jsonfile)

master_commands = []
saved_file = 'dummy.root'
files_to_ship = ['modifyCfg.py']

for i, seq in enumerate(sequence):
    release = seq['release'] 
    JOB_SCRIPT = SINGLE_JOB

    if release == 'local':
        scram_arch = os.environ['SCRAM_ARCH']
    else:
        scram_arch = seq['SCRAM_ARCH']

    run_script = seq['cfg']
    files_to_ship.append(str(run_script))
    files_to_ship.extend([str(x) for x in seq['local_files']])
    files_to_ship.extend([str(x) for x in seq['cmssw_files']])

    if args.crab:
        seq['local_files'] = [os.path.basename(x) for x in seq['local_files']]
        seq['cmssw_files'] = [os.path.basename(x) for x in seq['cmssw_files']]

    copy_arg = '\n'.join(['cp ${RUN_DIR}/%s ${CHAIN_DIR}/%s' % (x, os.path.basename(x)) for x in seq['local_files']])
    copy_arg += '\n'
    copy_arg += '\n'.join(['cp ${RUN_DIR}/%s ${CHAIN_DIR}/${CMSSW_RELEASE}/src/%s' % (x, os.path.basename(x)) for x in seq['cmssw_files']])

    events = -1
    inputFileArg = '--inputFile=step_%i.root' % (i - 1)
    outputFileArg = '--outputFile=step_%i.root' % (i)
    lumiOffsetArg = ''
    if i == 0:
        events = args.events
        inputFileArg = ''
        lumiOffsetArg = '--setLumiOffsets 500'
    if args.crab:
        run_script = os.path.basename(run_script)
    inputScript = os.path.join('${RUN_DIR}', run_script)
    outputScript = os.path.join('${CHAIN_DIR}', 'step_%i_cfg.py' % i)
    if 'output_module' in seq:
        outputModuleArg = '--outputModule %s' % seq['output_module']
    else:
        outputModuleArg = ''
    mod_arg = 'python ${RUN_DIR}/modifyCfg.py %s %s --events=%i --randomSeeds=${SEED} %s %s %s %s' % (inputScript, outputScript, events, inputFileArg, outputFileArg, outputModuleArg, lumiOffsetArg)

    return_files_arg = ''
    
    jobreport_name = 'JobReport_%i.xml' % i
    if 'save' in seq and seq['save'] == True:
        jobreport_name = 'FrameworkJobReport.xml'
        if 'output_files' not in seq:
            seq['output_files'] = []
        seq['output_files'].append(jobreport_name)
        seq['output_files'].append('step_%i.root' % i)
        saved_file = 'step_%i.root' % i
    if 'output_files' in seq:
        return_files_arg = '\n'.join(['cp ${CHAIN_DIR}/%s ${RUN_DIR}/%s' % (x, x) for x in set(seq['output_files'])])
    
    gen_script = JOB_SCRIPT % ({
        'SCRAM_ARCH': scram_arch,
        'CMSSW_RELEASE': release,
        'RUN_SCRIPT': outputScript,
        'FILE_COPY': copy_arg,
        'MOD_SCRIPT': mod_arg,
        'RETURN_FILES': return_files_arg,
        'JOB_REPORT': jobreport_name 
        })
    #print gen_script
    script_name = 'chain_step_%i_%s.sh' % (i, args.label)
    files_to_ship.append(script_name)

    with open(script_name, "w") as text_file:
        text_file.write(gen_script)

    master_commands.append('bash %s ${SEED}' % script_name)

master_script = MASTER_JOB % ({
    'CONTENT': '\n'.join(master_commands)
  })

master_script_name = 'run_chain_%s.sh' % args.label
with open(master_script_name, "w") as text_file:
    text_file.write(master_script)
os.system('chmod +x %s ' % (master_script_name))
files_to_ship.append(master_script_name)

cfg_other_settings = ''
real_cfg_name ='do_nothing_%s_cfg.py' % args.label

if args.crab and (args.crab_cfg is not None):
    d = {}
    execfile(args.crab_cfg, d)
    config = d['config']

    config.JobType.outputFiles = [saved_file]
    #config.JobType.numCores = 4
    #config.JobType.maxMemoryMB = 6000
    config.JobType.scriptExe = master_script_name
    config.JobType.psetName = real_cfg_name
    config.JobType.inputFiles = files_to_ship
    config.Data.unitsPerJob = args.events

    # Figure out if numCores has been set, and if so make sure
    # the same value is set in do_nothing_X_cfg.py,
    # otherwise crab will complain and refuse to submit
    if hasattr(config.JobType, 'numCores'):
        numCores = config.JobType.numCores
        cfg_other_settings = 'process.options.numberOfThreads=cms.untracked.uint32(%i)\nprocess.options.numberOfStreams=cms.untracked.uint32(0)' % numCores

    print config.__str__()
    with open('crab_generated_%s_cfg.py' % args.label, "w") as ofile:
        ofile.write(config.__str__())
        ofile.close()
    #os.system('chmod +x %s; ./%s' % (script_name, script_name))

if args.crab:
    real_cfg = DO_NOTHING_CFG % ({
        'FILENAME': saved_file,
        'OTHER_SETTINGS': cfg_other_settings
        })
    with open(real_cfg_name, "w") as ofile:
        ofile.write(real_cfg)
        ofile.close()
