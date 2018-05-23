import FWCore.ParameterSet.Config as cms
process = cms.Process("MAIN")
# import sys

# ################################################################
# # Read Options
# ################################################################
import FWCore.ParameterSet.VarParsing as parser
opts = parser.VarParsing ('analysis')
# #opts.register('file', 'root://xrootd.unl.edu//store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/0863B733-1A39-E611-AF47-0025905C53D8.root', parser.VarParsing.multiplicity.singleton,
# #opts.register('file', 'root://xrootd.unl.edu//store/mc/RunIISummer16MiniAODv2/SUSYGluGluToBBHToTauTau_M-1000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/4C466283-6BC0-E611-B3AE-001517FB25E4.root', parser.VarParsing.multiplicity.singleton,
# opts.register('file', 'root://xrootd.unl.edu//store/mc/RunIISummer16MiniAODv2/SUSYGluGluToBBHToTauTau_M-1000_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/130000/10B3D2AA-286C-E711-B57F-141877410B85.root', parser.VarParsing.multiplicity.singleton,               
# #opts.register('file', 'root://xrootd.unl.edu//store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/000C6E52-8BEC-E611-B3FF-0025905C42FE.root',parser.VarParsing.multiplicity.singleton,
# #opts.register('file', 'root://xrootd.unl.edu//store/data/Run2016F/SingleMuon/MINIAOD/PromptReco-v1/000/277/932/00000/084865EB-1859-E611-BDA7-02163E011A89.root', parser.VarParsing.multiplicity.singleton,
# #opts.register('file', 'root://xrootd.unl.edu//store/data/Run2016H/SingleMuon/MINIAOD/PromptReco-v2/000/281/265/00000/28861171-6E82-E611-9CAF-02163E0141FA.root', parser.VarParsing.multiplicity.singleton,
# #opts.register('file', 'root://xrootd.unl.edu//store/data/Run2016H/SingleElectron/MINIAOD/PromptReco-v3/000/284/036/00000/1CBE1DEB-589F-E611-ABBB-02163E0143B5.root', parser.VarParsing.multiplicity.singleton,
# parser.VarParsing.varType.string, "input file")
# opts.register('globalTag', '80X_mcRun2_asymptotic_2016_TrancheIV_v7', parser.VarParsing.multiplicity.singleton,
# #opts.register('globalTag', '80X_dataRun2_2016SeptRepro_v7', parser.VarParsing.multiplicity.singleton,
#     parser.VarParsing.varType.string, "global tag")
# opts.register('isData', 0, parser.VarParsing.multiplicity.singleton,
#     parser.VarParsing.varType.int, "Process as data?")

opts.register('cores', 1, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Number of cores/threads")
opts.register('input', 'root://xrootd.unl.edu//store/data/Run2016H/Tau/MINIAOD/PromptReco-v3/000/284/036/00000/36B9BD65-5B9F-E611-820B-02163E0126D3.root', parser.VarParsing.multiplicity.singleton, parser.VarParsing.varType.string, "input file")


# #opts.register('release', '7412MINIAOD', parser.VarParsing.multiplicity.singleton,
# opts.register('release', '80XMINIAOD', parser.VarParsing.multiplicity.singleton,
#     parser.VarParsing.varType.string, "Release label")
# opts.register('doHT', 0, parser.VarParsing.multiplicity.singleton,
#     parser.VarParsing.varType.int, "Store HT and number of outgoing partons?")
# opts.register('isReHLT', 1, parser.VarParsing.multiplicity.singleton,
#     parser.VarParsing.varType.int, "Process as reHLT sample?")
# opts.register('LHEWeights', False, parser.VarParsing.multiplicity.singleton,
#     parser.VarParsing.varType.bool, "Produce LHE weights for sample")
# opts.register('LHETag', 'externalLHEProducer', parser.VarParsing.multiplicity.singleton,
#     parser.VarParsing.varType.string, "Input tag for LHE weights")



opts.parseArguments()
# infile      = opts.file
# if not infile: infile = "file:/tmp/file.root"
# isData      = opts.isData
# tag         = opts.globalTag
# release     = opts.release
# doLHEWeights = opts.LHEWeights
# if not isData:
#   doHT     = opts.doHT
#   isReHLT  = opts.isReHLT
# else:
#   doHT     = 0
#   isReHLT  = 0
# #isEmbedded  = opts.isEmbedded
# #isTandP     = opts.isTandP
# #isZStudy    = opts.isZStudy
# #isPhys14    = opts.isPhys14

# if not release in ["76X", "80XMINIAOD"]:
#   print 'Release not recognised, exiting!'
#   sys.exit(1)
# print 'release     : '+release
# print 'isData      : '+str(isData)
# print 'globalTag   : '+str(tag)
# print 'doHT        : '+str(doHT)
# print 'isReHLT     : '+str(isReHLT)

################################################################
# Standard setup
################################################################
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

################################################################
# Message Logging, summary, and number of events
################################################################
isMC = True

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(1000)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    wantSummary=cms.untracked.bool(True)
)

################################################################
# Input files and global tags
################################################################
process.load("CondCore.CondDB.CondDB_cfi")
from CondCore.CondDB.CondDB_cfi import *

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
    # /store/mc/RunIISummer16MiniAODv2/WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/FEB2F873-C1D8-E611-8AAC-02163E012A61.root 
    # 'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/94X_mc2017_realistic_v10_ext1-v1/00000/0000BD66-99F4-E711-97DF-24BE05C33C22.root',
    # 'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAOD/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/001F8E98-4D05-E811-91A5-02163E019B9C.root'
    # '/store/mc/RunIIFall17MiniAODv2/Z1JetsToNuNu_M-50_LHEZpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/00430106-6B42-E811-BA09-001F29089F7E.root'
   opts.input 
))

# 80X_mcRun2_asymptotic_2016_TrancheIV_v8
# 94X_mc2017_realistic_v10
process.GlobalTag.globaltag = cms.string('80X_mcRun2_asymptotic_2016_TrancheIV_v8')


process.options.numberOfThreads = cms.untracked.uint32(opts.cores)
process.options.numberOfStreams = cms.untracked.uint32(opts.cores)


process.acMuonProducer = cms.EDProducer('AcornMuonProducer',
    input=cms.InputTag("slimmedMuons"),
    branch=cms.string('muons'),
    inputVertices=cms.InputTag('offlineSlimmedPrimaryVertices'),
    select=cms.vstring('keep .* p4=12 pfIso.*=12')
)


#### Adding the photon ID
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process, idmod, setupVIDPhotonSelection)


process.acPhotonProducer = cms.EDProducer('AcornPhotonProducer',
    input=cms.InputTag("slimmedPhotons"),
    branch=cms.string('photons'),
    select=cms.vstring('keep .* p4=12'),
    phoLooseIdMap=cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-loose"),
    phoMediumIdMap=cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium"),
    phoTightIdMap=cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-tight")
    # phoMediumIdFullInfoMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium")
)

process.acGenParticleProducer = cms.EDProducer('AcornGenParticleProducer',
    input=cms.InputTag("prunedGenParticles"),
    branch=cms.string('genParticles'),
    select=cms.vstring('keep .* p4=12', 'drop mothers')
)

process.acLHEParticleProducer = cms.EDProducer('AcornLHEParticleProducer',
    input=cms.InputTag("externalLHEProducer"),
    branch=cms.string('lheParticles'),
    select=cms.vstring('keep .* p4=12')
)

process.acPileupInfoProducer = cms.EDProducer('AcornPileupInfoProducer',
    input=cms.InputTag("slimmedAddPileupInfo"),
    branch=cms.string('pileupInfo'),
    select=cms.vstring('keep .*')
)

process.acEventInfoProducer = cms.EDProducer('AcornEventInfoProducer',
    lheProducer=cms.InputTag("externalLHEProducer"),
    generator=cms.InputTag("generator"),
    includeLHEWeights=cms.bool(isMC),
    includeGenWeights=cms.bool(isMC),
    branch=cms.string('eventInfo'),
    select=cms.vstring(
        'keep .*',
        'drop lheweights:.*',
        'keep lheweights:(renscfact|facscfact|muR|muF|mur|muf|MUR|MUF).*=10',
        'keep lheweights:lhapdf.306[0-9][0-9][0-9]=10',
        'keep lheweights:PDF.306000=10',
        'keep lheweights:dim6=10',
        'keep lheweights:NNPDF31_nnlo_hessian_pdfas=10')
)

process.acEventProducer = cms.EDProducer('AcornEventProducer')

process.p = cms.Path(
    process.egmPhotonIDSequence +
    process.acMuonProducer +
    process.acPhotonProducer +
    process.acLHEParticleProducer +
    process.acGenParticleProducer +
    process.acPileupInfoProducer +
    process.acEventInfoProducer +
    process.acEventProducer)

# process.schedule = cms.Schedule(process.patTriggerPath, process.p)
process.schedule = cms.Schedule(process.p)
# print process.dumpPython()
