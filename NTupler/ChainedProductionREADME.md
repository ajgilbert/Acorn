# Chained production

The script `cmsRunChain.py` is used to configure and generate a set of shell scripts that will take care of running each stage of the job (setting up the appropriate CMSSW area, modifying the cmsRun PSets, etc.). The configuration is stored in a simple json file. Typical usage:

    python cmsRunChain.py vhbb_private.json --events 2000 [--crab]

By default the scripts produced can be used for testing locally. Add the option `--crab` to generate versions that will work in a crab job. The json config file looks like this:

```json
[
  {
    "release": "CMSSW_9_3_6_patch1",
    "SCRAM_ARCH": "slc6_amd64_gcc630",
    "cfg": "vhbb-miniaod-2017/HIG-RunIIFall17wmLHEGS-00620_1_cfg.py",
    "local_files": [],
    "cmssw_files": [],
    "output_module": "RAWSIMoutput"
  },
  {
    "release": "CMSSW_9_4_0_patch1",
    "SCRAM_ARCH": "slc6_amd64_gcc630",
    "cfg": "vhbb-miniaod-2017/HIG-RunIIFall17DRPremix-00742_1_cfg.py",
    "local_files": [],
    "cmssw_files": []
  },
  {
    "release": "CMSSW_9_4_0_patch1",
    "SCRAM_ARCH": "slc6_amd64_gcc630",
    "cfg": "vhbb-miniaod-2017/HIG-RunIIFall17DRPremix-00742_2_cfg.py",
    "local_files": [],
    "cmssw_files": []
  },
  {
    "release": "CMSSW_9_4_6_patch1",
    "SCRAM_ARCH": "slc6_amd64_gcc630",
    "cfg": "vhbb-miniaod-2017/HIG-RunIIFall17MiniAODv2-00605_1_cfg.py",
    "local_files": [],
    "cmssw_files": [],
    "save": true,
    "output_files": ["step_3.root", "FrameworkJobReport.xml"]
  }
]
```

Each block gives the basic information needed to run the job - the CMSSW release and SCRAM_ARCH to use, and the path to the cmsRun config. The other options are:

 * `local_files`: a list of files that should be copied into the area where cmsRun will be executed (NB: this will not be inside the CMSSW directory that is created for the job, see the next option)
 * `cmssw_files`: a list of files that should be copied into the `$CMSSW_BASE/src` directory at the start of the job
 * `output_module`: when each step runs it will automatically modify the cmsRun PSet to make sure the output file of one step is consistent with the input file of the next. It does this by searching for an OutputModule in the config. In some cases more than one can be found, so here we specify explicitly which one we want to use with the subsequent job.
 * `save`: The output from this step should be published by crab - NB can do this for one step only.
 * `output_files`: List of files produced by this step that should be returned to where master job is running. For the step where the output will be published this should include the ROOT output file and the `FrameworkJobReport.xml` file that will be created by cmsRun.

Note that the input cmsRun configs do have to be modified slightly compared to what is retrieved from McM:

 * A placeholder (`#{EDITS}`) for injecting additional parameter settings should be added just after the end of the customize part:

```py
#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions
#{EDITS}

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
```

 * The multithreading settings should be hardcoded here. For many McM workflows they may have been set already, for some they might not have been set at all. E.g. to use four cores add:

```py
process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(0)
```

The output from running `cmsRunChain.py` is the main job script, `run_chain.sh`, and a separate script for each step i, `chain_step_i.sh`. To test locally can simply do:

    bash run_chain.sh <RNG seed>


## Running with crab

An example crab config (`crab_vhbb.py`) is included. There are currently a few considerations

 * The list of additional input files has to be set manually, and should include the master job script (`run_chain.sh`), the scripts for each step, the `modifyCfg.py` script, the input cmsRun configs, and any additional files that were specified in the json config above.
 * The `unitsPerJob` and `totalUnits` will be used by crab as normal to decide the number of jobs to submit, but note that the number of events per job is always fixed by the option to `cmsRunChain.py` above.
 * The `numCores` should be set appropriately, and also needs to be set to the same value in the dummy `do_nothing_cfg.py` config, because crab will not allow a mismatch here.
