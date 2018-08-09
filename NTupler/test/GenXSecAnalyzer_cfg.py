# Auto generated configuration file
# using: 
# Revision: 1.381.2.7 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/EightTeV/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_tauola_cff.py --step GEN --beamspot Realistic8TeVCollision --conditions START52_V9::All --pileup NoPileUp --datamix NODATAMIXER --eventcontent RAWSIM --datatier GEN -n -1 --python_filename=Temp_Hadronizer_5498_1.py --filein root://eoscms//eos/cms//store/lhe/5498/DY1JetsToLL_M-50_8TeV-madgraph_10001.lhe --no_output --no_exec
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register('events', -1, VarParsing.multiplicity.singleton,
    VarParsing.varType.int, "Number of events")
options.register('input', 'root://xrootd.unl.edu//store/data/Run2016H/Tau/MINIAOD/PromptReco-v3/000/284/036/00000/36B9BD65-5B9F-E611-820B-02163E0126D3.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, "input file")
options.parseArguments()

process = cms.Process('ANA')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.events)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source(
    "PoolSource",
    fileNames  = cms.untracked.vstring(options.input.split(',')),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)


process.dummy2 = cms.EDAnalyzer("GenXSecAnalyzer")



# Path and EndPath definitions
process.ana = cms.Path(process.dummy2)
# Schedule definition
process.schedule = cms.Schedule(process.ana)
