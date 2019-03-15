import os
import json
import stat
import copy
import glob
import datetime
from math import ceil
from functools import partial
from multiprocessing import Pool

PROVISION_CMSSW_AREA = """CMSSW_RELEASE=%(CMSSW_RELEASE)s

if [ -d ${CMSSW_RELEASE} ]; then
  echo ">> Directory ${CMSSW_RELEASE} already exists, quiting"
  exit 1
fi

xrdcp %(REMOTE_TARBALL)s cmssw.tar.gz

echo ">> Setting up release area for ${CMSSW_RELEASE} and arch ${SCRAM_ARCH}"
if [ ! -d ${CMSSW_RELEASE} ]; then
  scram project CMSSW ${CMSSW_RELEASE}
fi

tar -xf cmssw.tar.gz -C ${CMSSW_RELEASE}/
rm cmssw.tar.gz
"""

JOB_PREFIX = """#!/bin/sh
export INITIALDIR=${PWD}
set -o pipefail
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=%(SCRAM_ARCH)s

%(PROVISION_CMSSW_AREA)s
cd %(CMSSW_BASE)s/src
eval `scramv1 runtime -sh`
cd %(PWD)s
"""

JOB_NAF_PREFIX = """#!/bin/sh
export INITIALDIR=${PWD}
set -o pipefail
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd %(CMSSW_BASE)s/src
export SCRAM_ARCH=%(SCRAM_ARCH)s
eval `scramv1 runtime -sh`
cd %(PWD)s
"""

CONDOR_TEMPLATE = """executable = %(EXE)s
arguments = $(ProcId)
output                = %(TASKDIR)s%(TASK)s.$(ClusterId).$(ProcId).out
error                 = %(TASKDIR)s%(TASK)s.$(ClusterId).$(ProcId).err
log                   = %(TASKDIR)s%(TASK)s.$(ClusterId).log

# Send the job to Held state on failure.
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)

# Periodically retry the jobs every 10 minutes, up to a maximum of 5 retries.
# periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)

%(EXTRA)s
queue %(NUMBER)s

"""

ERROR_CAPTURE = """
function error_exit
{
    mv %s %s.$1
    exit $1
}

trap 'error_exit $?' ERR
"""


def run_command(dry_run, command):
    if not dry_run:
        print '>> ' + command
        return os.system(command)
    else:
        print '[DRY-RUN]: ' + command


class Jobs:
    description = 'Simple job submission system'

    def __init__(self):
        self.job_queue = []
        self.args = None
        self.job_mode = 'interactive'
        self.parallel = 1
        self.merge = 1
        self.task_name = 'combine_task'
        self.dry_run = False
        self.bopts = ''  # batch submission options
        self.tracking = False
        self.task_dir = ''

    def attach_job_args(self, group):
        group.add_argument('--job-mode', default=self.job_mode, choices=[
                           'interactive', 'script', 'lxbatch', 'condor', 'NAF', 'ts'], help='Task execution mode')
        group.add_argument('--task-name', default=self.task_name,
                           help='Task name, used for job script and log filenames for batch system tasks')
        group.add_argument('--dir', default=self.task_dir,
                           help='Area for creating the job scripts')
        group.add_argument('--parallel', type=int, default=self.parallel,
                           help='Number of jobs to run in parallel [only affects interactive job-mode]')
        group.add_argument('--merge', type=int, default=self.merge,
                           help='Number of jobs to run in a single script [only affects batch submission]')
        group.add_argument('--dry-run', action='store_true',
                           help='Print commands to the screen but do not run them')
        group.add_argument('--cwd', type=int, default=1,
                           help='Switch to the submission directory within the job')
        group.add_argument('--sub-opts', default=self.bopts,
                           help='Options for batch submission')
        group.add_argument('--memory', type=int,
                           help='Request memory for job [MB]')
        group.add_argument('--tracking', nargs='?', default=False, const='short',
                           help='Track job status (if applicable)')
        group.add_argument('--ship', type=str, default=None,
                           help='JSON config for shipping local CMSSW area to the remote hosts')
        group.add_argument('--ship-tarball', type=str, default=None,
                           help='Use this filename for the CMSSW tarball')

    def set_args(self, args):
        self.args = args
        self.job_mode = self.args.job_mode
        self.task_name = self.args.task_name
        self.parallel = self.args.parallel
        self.merge = self.args.merge
        self.dry_run = self.args.dry_run
        self.bopts = self.args.sub_opts
        self.memory = self.args.memory
        self.tracking = self.args.tracking
        self.task_dir = self.args.dir
        self.ship = self.args.ship
        self.ship_tarball = self.args.ship_tarball
        if self.ship is not None and self.ship_tarball is None:
            self.ship_tarball = '%s_%s.tar.gz' % (os.environ['CMSSW_VERSION'], datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
        self.provision_cmd = ''
        # if self.dry_run:
        #     self.tracking = False

    def file_len(self, fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
            return i + 1

    def add_filelist_split_jobs(self, prog, cfg, files_per_job, output_cfgs=['output'], tempdir=''):
        if 'filelist' in cfg:
            filelist = cfg['filelist']
            nfiles = self.file_len(filelist)
        elif 'filelists' in cfg:
            nfiles = 0
            for f in cfg['filelists']:
                nfiles += self.file_len(f)
        if files_per_job <= 0:
            njobs = 1
        else:
            njobs = int(ceil(float(nfiles)/float(files_per_job)))
        for n in xrange(njobs):
            final_cfg = copy.deepcopy(cfg)
            for item in output_cfgs:
                filename, extension = os.path.splitext(cfg[item])
                final_cfg[item] = filename + ('_%i' % n) + extension
            if tempdir:
                copy_cmds = []
                cp_cmd = 'cp'
                if final_cfg['outdir'].startswith('root://'):
                    cp_cmd = 'xrdcp --force'
                for item in output_cfgs:
                    copy_cmds.append("""%s %s/%s %s/%s""" % (cp_cmd, tempdir, final_cfg[item], final_cfg['outdir'], final_cfg[item]))
                final_cfg['outdir']=tempdir
            final_cfg['file_offset'] = n
            final_cfg['file_step'] = njobs
            cmd = """%s '%s'""" % (prog, json.dumps(final_cfg))
            if tempdir:
                cmd = """%s '%s'; %s""" % (prog, json.dumps(final_cfg), '; '.join(copy_cmds))
            self.job_queue.append(cmd)
            # print cmd

    def create_job_script(self, commands, script_filename, do_log=False):
        fname = script_filename
        logname = script_filename.replace('.sh', '.log')
        DO_JOB_PREFIX = JOB_PREFIX if not self.job_mode == 'NAF' else JOB_NAF_PREFIX
        DO_JOB_PREFIX = DO_JOB_PREFIX % ({
          'CMSSW_BASE': os.environ['CMSSW_BASE'],
          'SCRAM_ARCH': os.environ['SCRAM_ARCH'],
          'PWD': (os.environ['PWD'] if self.args.cwd else '${INITIALDIR}'),
          'PROVISION_CMSSW_AREA': ''
        })

        with open(fname, "w") as text_file:
            text_file.write(DO_JOB_PREFIX)
            if self.tracking:
                full_path = os.path.abspath(script_filename)
                text_file.write('mv %s %s\n' % (
                        full_path.replace('.sh', '.status.submitted'),
                        full_path.replace('.sh', '.status.running')
                    ))
                text_file.write(ERROR_CAPTURE % (
                        full_path.replace('.sh', '.status.running'),
                        full_path.replace('.sh', '.status.error')
                    ))
            for i, command in enumerate(commands):
                tee = 'tee' if i == 0 else 'tee -a'
                log_part = '\n'
                if do_log: log_part = ' 2>&1 | %s ' % tee + logname + log_part
                text_file.write(command + log_part)
            if self.tracking:
                full_path = os.path.abspath(script_filename)
                text_file.write('mv %s %s\n' % (
                        full_path.replace('.sh', '.status.running'),
                        full_path.replace('.sh', '.status.done')
                    ))
        st = os.stat(fname)
        os.chmod(fname, st.st_mode | stat.S_IEXEC)
        print 'Created job script: %s' % script_filename
        if self.tracking and not self.dry_run:
            trackname = script_filename.replace('.sh', '.status.created')
            open(trackname, 'a').close()

    def read_job_status(self, script_filename):
        full_path = os.path.abspath(script_filename)
        status_files = glob.glob(full_path.replace('.sh', '.status.*'))
        statuses = []
        for status_f in status_files:
            pos = status_f.rfind('status')
            status = (status_f[pos:]).split('.')[1:]
            statuses.append(status)
        # print statuses
        return statuses

    def flush_queue(self):
        if self.job_mode == 'interactive':
            pool = Pool(processes=self.parallel)
            result = pool.map(
                partial(run_command, self.dry_run), self.job_queue)
        # flush_queue might be called more than once but only need
        # to upload the CMSSW tarball once. After the first call
        # self.provision_cmd will no longer be an empty string, so
        # use this to check.
        if self.ship is not None and self.provision_cmd == '':
            self.upload_tarball()
        script_list = []
        status_result = {}
        njobs = 0
        if self.job_mode in ['script', 'lxbatch', 'NAF', 'ts']:
            for i, j in enumerate(range(0, len(self.job_queue), self.merge)):
                njobs += 1
                script_name = 'job_%s_%i.sh' % (self.task_name, i)
                if self.task_dir is not '':
                    script_name = os.path.join(self.task_dir, script_name)
                status = 'unknown'
                if self.tracking:
                    statuses = self.read_job_status(script_name)
                    if len(statuses) == 0:
                        status = 'new'
                    elif len(statuses) == 1:
                        status = statuses[0][0]
                    else:
                        status = 'confused'
                    if status not in status_result:
                       status_result[status] = []
                    status_result[status].append(script_name)
                    if len(statuses) > 0:
                        # print '%s appears to already be in progress, skipping...' % script_name
                        continue

                # each job is given a slice from the list of combine commands of length 'merge'
                # we also keep track of the files that were created in case submission to a
                # batch system was also requested
                self.create_job_script(
                    self.job_queue[j:j + self.merge], script_name, self.job_mode in ['script', 'ts'])
                script_list.append(script_name)
        if self.job_mode == 'lxbatch':
            for script in script_list:
                full_script = os.path.abspath(script)
                logname = full_script.replace('.sh', '_%J.log')
                if self.tracking and not self.dry_run:
                    os.rename(full_script.replace('.sh', '.status.created'), full_script.replace('.sh', '.status.submitted'))
                run_command(self.dry_run, 'bsub -o %s %s %s' % (logname, self.bopts, full_script))
        if self.job_mode == 'NAF':
            for script in script_list:
                full_script = os.path.abspath(script)
                logname = full_script.replace('.sh', '_$JOB_ID.log')
                if self.tracking and not self.dry_run:
                    os.rename(full_script.replace('.sh', '.status.created'), full_script.replace('.sh', '.status.submitted'))
                run_command(self.dry_run, """qsub -j y -o '%s' %s %s""" % (logname, self.bopts, full_script))
        if self.job_mode == 'ts':
            for script in script_list:
                full_script = os.path.abspath(script)
                if self.tracking and not self.dry_run:
                    os.rename(full_script.replace('.sh', '.status.created'), full_script.replace('.sh', '.status.submitted'))
                run_command(self.dry_run, 'ts bash -c "eval %s"' % (full_script))
        if self.job_mode == 'condor':
            outscriptname = 'condor_%s.sh' % self.task_name
            subfilename = 'condor_%s.sub' % self.task_name
            if self.task_dir is not '':
                outscriptname = os.path.join(self.task_dir, outscriptname)
                subfilename = os.path.join(self.task_dir, subfilename)
            print '>> condor job script will be %s' % outscriptname
            outscript = open(outscriptname, "w")
            job_run_dir = os.environ['PWD'] if self.args.cwd else '${INITIALDIR}'
            if self.ship:
                job_run_dir = '%s/%s/%s' % ('${INITIALDIR}', '${CMSSW_RELEASE}', os.path.relpath(os.environ['PWD'], os.environ['CMSSW_BASE']))
            DO_JOB_PREFIX = JOB_PREFIX % ({
              'CMSSW_BASE': os.environ['CMSSW_BASE'] if not self.ship else '${CMSSW_RELEASE}',
              'SCRAM_ARCH': os.environ['SCRAM_ARCH'],
              'PWD': job_run_dir,
              'PROVISION_CMSSW_AREA': self.provision_cmd
            })
            outscript.write(DO_JOB_PREFIX)
            jobs = 0
            for i, j in enumerate(range(0, len(self.job_queue), self.merge)):
                outscript.write('\nif [ $1 -eq %i ]; then\n' % jobs)
                jobs += 1
                for line in self.job_queue[j:j + self.merge]:
                    newline = line
                    outscript.write('  ' + newline + '\n')
                outscript.write('fi')
            outscript.close()
            st = os.stat(outscriptname)
            os.chmod(outscriptname, st.st_mode | stat.S_IEXEC |stat.S_IXGRP | stat.S_IXOTH)
            subfile = open(subfilename, "w")
            condor_settings = CONDOR_TEMPLATE % {
              'EXE': outscriptname,
              'TASK': self.task_name,
              'TASKDIR': os.path.join(self.task_dir, ''),
              'EXTRA': self.bopts.decode('string_escape'),
              'NUMBER': jobs
            }
            print '>> Created %i jobs to process %i tasks' % (jobs, len(self.job_queue))
            subfile.write(condor_settings)
            subfile.close()
            run_command(self.dry_run, 'condor_submit %s' % (subfilename))

        if self.job_mode in ['lxbatch', 'NAF', 'ts'] and self.tracking:
            print '>> Status summary: %s' % self.task_name
            for status in status_result:
                counter = '[%i/%i]' % (len(status_result[status]), njobs)
                print '%20s %10s' % (status, counter)
                if self.tracking == 'long':
                    for f in status_result[status]:
                        print ' '*20 + '%s' % f

        del self.job_queue[:]

    def upload_tarball(self):
        # Read the ship JSON file
        with open(self.ship) as jsonfile:
            ship_cfg = json.load(jsonfile)
        remote_dir = ship_cfg['remote_dir']
        manifest = ship_cfg['manifest']
        run_command(self.dry_run, 'pushd %(CMSSW_BASE)s; tar --exclude=".*" -zcf %(TARBALL)s %(MANIFEST)s; xrdcp --force %(TARBALL)s %(REMOTEDIR)s/%(TARBALL)s; rm %(TARBALL)s; popd' % {
            'CMSSW_BASE': os.environ['CMSSW_BASE'],
            'TARBALL': self.ship_tarball,
            'MANIFEST': ' '.join(manifest),
            'REMOTEDIR': remote_dir
            })
        self.provision_cmd = PROVISION_CMSSW_AREA % {
            'CMSSW_RELEASE': os.environ['CMSSW_VERSION'],
            'REMOTE_TARBALL': '%s/%s' % (remote_dir, self.ship_tarball)
        }
