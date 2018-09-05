# Acorn Framework

This page will document the main tools available in the framework, in roughly the order they would be used in an analysis workflow.

## NTuple production

### CMSSW modules and config

An example cmsRun config can be found in `NTupler/test/wgamma_cfg.py`.


Most object collection plugins are dervied from a common EDProducer: `AcornBaseProducer`. This is a template class, where the template parameter should be the type of collection our derived producer will make:

```c++
        class AcornMuonProducer : public AcornBaseProducer<std::vector<ac::Muon>>
```

The `AcornBaseProducer` handles a few common tasks automatically:

 - Provides a `branch` setting in the PSet to specify the branch name in the output TTree, e.g. `branch=cms.string('photons')`.
 - Takes care of creating a pointer to the output collection, and connecting this to the output TTree (which is itself created by the `AcornEventProducer`).
 - Provides a `select` setting that can be used to specify exactly which variables should be saved in the output object, and optionally allows the numerical precision to be truncated (inspired by miniAOD/nanoAOD). The syntax is similar to that used in the GenParticlePruner: a series of `keep` and `drop` rules are given in order, which is then checked for each variable. The basic pattern is `keep/drop regex[=N]`, where the `=N` is optional, and specifies the number of bits on the mantissa to keep for floating point variables.
 - Works correctly as part of the CMSSW multi-threaded framework (if enabled in the config via `process.options.numberOfThreads = cms.untracked.uint32(N)` and `process.options.numberOfStreams = cms.untracked.uint32(N)`)

### Configure a large-scale production

Can use a json file to keep track of all the samples to be processed, and what cmsRun options should be used with which sample. See `NTupler/test/samples/wgamma_2016.json` for an example. The basic format is:

```json
{
  "samples": {
    "SingleMuon-2016G": {
      "dataset": "/SingleMuon/Run2016G-03Feb2017-v1/MINIAOD",
      "attributes": ["submit"],
      "events": 149916849,
      "config": "data_reproc"
    },
    "SingleMuon-2016H-ver2": {
      "dataset": "/SingleMuon/Run2016H-03Feb2017_ver2-v1/MINIAOD",
      "attributes": ["submit"],
      "events": 169642135,
      "config": "data_prompt"
    }
  },
  "configs": {
    "data_reproc": ["year=2016", "isData=1", "globalTag=80X_dataRun2_2016SeptRepro_v7", "cores=1"],
    "data_prompt": ["year=2016", "isData=1", "globalTag=80X_dataRun2_Prompt_v16", "cores=1"]
  }
}
```
Each sample gives the dataset name and a config label, which is mapped to the set of options at the bottom. It is not essential to fill the number of events, but could be useful for book-keeping later. The list of attributes are to be used by other scripts which might want to select some subset of the entries in the json to act on.

#### cmsRunTest.py: Testing the cmsRun config on a sample

Example usage:

```bash
python test/cmsRunTest.py test/samples/wgamma_2016.json --samples ZZTo4L-amcatnloFXFX test/wgamma_cfg.py
```

The first option is the chosen json file, the `--samples` option is a comma separated list of samples to check from that file. All remaining arguments will be passed through to `cmsRun` (here we just give the cmsRun config). This script will query DBS to find a file from the sample, and will automatically add an option `input=[filename]` to the cmsRun config - so important to make sure your config will accept this option.

#### Running the GenXSecAnalyser

A barebones cmsRun config that just runs the GenXSecAnalyzer is available in `test/GenXSecAnalyzer_cfg.py`. Can be used in conjunction with `cmsRunTest.py`:

```bash
python test/cmsRunTest.py test/samples/wgamma_2016.json --samples WZTo1L3Nu-amcatnloFXFX --no-cfg --nfiles 3 test/GenXSecAnalyzer_cfg.py maxEvents=-1
```

Here the `--no-cfg` option is added to suppress the automatic addition of the cmsRun command line options from the json file, and `--nfiles X` lets us process more files from the sample to increase the stats for the cross section info.

#### crabSubmit.py: Submitting with crab

A script for submitting jobs with crab:

```bash
python test/crabSubmit.py test/samples/wgamma_2016.json -p test/wgamma_cfg.py -l wgamma_2016_v1 -u 250000 --attribute submitgenonly -v 1 [--submit]
```

Where the first argument is the json file, `-p` is the cmsRun config, `-l` is a label for the production (`config.General.workArea` and `config.Data.outLFNDirBase` will use this), `-u` sets the unitsPerJob (i.e. number of events per job) and `-v` is the verbosity of the script. The `--attribute` option can be used to submit only the samples that have that particular attribute assigned in the json file.

**FIXME**: the output path is current hard-coded - needs to be an option!

#### crabFilelist.py: Making a filelist

The script only needs the crab directory, from this it will check all expected output files are present and if so will produce the filelist:

```bash
python ../test/crabFilelist.py crab_wgamma_2016_v1_WGToLNuG-EFT-madgraphMLM
```

Additional options can be specified: `--replace-xrootd` can be used to set the local xrootd server, otherwise a global redirector will be used, making access slightly slower. The default is to replace with the CERN EOS server `eoscms.cern.ch`, so should override this at other sites. The `--prefix` option adds the given string to the beginning of the filelist (i.e. `[prefix]wgamma_2016_v1_WGToLNuG-EFT-madgraphMLM.txt` for the above example).


