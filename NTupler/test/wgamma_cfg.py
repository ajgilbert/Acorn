import FWCore.ParameterSet.Config as cms
process = cms.Process("MAIN")
# import sys

# ################################################################
# # Read Options
# ################################################################
import FWCore.ParameterSet.VarParsing as parser
opts = parser.VarParsing ('analysis')
opts.register('globalTag', '80X_dataRun2_2016SeptRepro_v7', parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.string, "global tag")
opts.register('events', 1000, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Number of events")
opts.register('isData', 0, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Process as data?")
opts.register('cores', 1, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Number of cores/threads")
opts.register('input', 'root://xrootd.unl.edu//store/data/Run2016H/Tau/MINIAOD/PromptReco-v3/000/284/036/00000/36B9BD65-5B9F-E611-820B-02163E0126D3.root', parser.VarParsing.multiplicity.singleton, parser.VarParsing.varType.string, "input file")
opts.register('year', '2016', parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.string, "Year label")

opts.parseArguments()
isData = bool(opts.isData)
isMC = not isData

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
process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(opts.events)
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

# Example files:
# 2016 MC: /store/mc/RunIISummer16MiniAODv2/WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/FEB2F873-C1D8-E611-8AAC-02163E012A61.root
# 2017 MC: /store/mc/RunIIFall17MiniAODv2/Z1JetsToNuNu_M-50_LHEZpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/00430106-6B42-E811-BA09-001F29089F7E.root
# 2018 MC: /store/mc/RunIISpring18MiniAOD/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/100X_upgrade2018_realistic_v10_ext1-v1/30000/00F3DB63-1D25-E811-B003-0025901D08B2.root
# 2017 data: /store/data/Run2017E/SingleMuon/MINIAOD/31Mar2018-v1/00000/000D53C5-9D39-E811-A39C-0025905B85A0.root
# 2018 data: /store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v2/000/316/239/00000/06C61F62-3759-E811-A213-02163E017F4E.root
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
   opts.input
))

process.GlobalTag.globaltag = cms.string(opts.globalTag)


process.options.numberOfThreads = cms.untracked.uint32(opts.cores)
process.options.numberOfStreams = cms.untracked.uint32(opts.cores)


process.selectedMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string("pt > 20 & abs(eta) < 2.6")
)

process.selectedPhotons = cms.EDFilter("PATPhotonRefSelector",
    src = cms.InputTag("slimmedPhotons"),
    cut = cms.string("pt > 20 & abs(eta) < 3.5")
)

process.acMuonProducer = cms.EDProducer('AcornMuonProducer',
    input=cms.InputTag("selectedMuons"),
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
if opts.year == '2016':
    my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff']
    photon_loose_id = "egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-loose"
    photon_medium_id = "egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium"
    photon_tight_id = "egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-tight"
else:
    my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V1_TrueVtx_cff']
    photon_loose_id = 'egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V1-loose'
    photon_medium_id = 'egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V1-medium'
    photon_tight_id = 'egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V1-tight'
#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process, idmod, setupVIDPhotonSelection)

process.acPhotonProducer = cms.EDProducer('AcornPhotonProducer',
    input=cms.InputTag("selectedPhotons"),
    branch=cms.string('photons'),
    select=cms.vstring('keep .* p4=12'),
    phoLooseIdMap=cms.InputTag(photon_loose_id),
    phoMediumIdMap=cms.InputTag(photon_medium_id),
    phoTightIdMap=cms.InputTag(photon_tight_id)
    # phoMediumIdFullInfoMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium")
)


process.acMCSequence = cms.Sequence(
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.selectedGenParticles = cms.EDProducer("GenParticlePruner",
    src=cms.InputTag("prunedGenParticles"),
    select=cms.vstring(
        "drop  *",
        "keep isPromptFinalState()",
        "keep++ abs(pdgId) == 15 && statusFlags().isPrompt()"
    )
)

process.acGenParticleProducer = cms.EDProducer('AcornGenParticleProducer',
    input=cms.InputTag("selectedGenParticles"),
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
    select=cms.vstring('keep .*'),
    minBx=cms.int32(0),
    maxBx=cms.int32(0)
)

if isMC:
    process.acMCSequence += cms.Sequence(
        process.selectedGenParticles +
        process.acGenParticleProducer +
        process.acLHEParticleProducer +
        process.acPileupInfoProducer
    )

hlt_paths = [
    # 'HLT_IsoMu22_v',
    # 'HLT_IsoTkMu22_v',
    'HLT_IsoMu24_v',
    'HLT_IsoTkMu24_v',
    'HLT_IsoMu27_v'
]

process.acTriggerObjectSequence = cms.Sequence(
    # process.icTriggerPathProducer
)

if opts.year == '2016':
    trigger_objects = 'selectedPatTrigger'
else:
    trigger_objects = 'slimmedPatTrigger'

for path in hlt_paths:
    shortname = path[4:-2]  # drop the HLT_ and _v parts
    setattr(process, 'ac_%s_ObjectProducer' % shortname, cms.EDProducer('AcornTriggerObjectProducer',
        input=cms.InputTag(trigger_objects),
        triggerResults=cms.InputTag('TriggerResults', '', 'HLT'),
        hltConfigProcess=cms.string('HLT'),
        branch=cms.string('triggerObjects_%s' % shortname),
        hltPath= cms.string(path),
        storeIfFired=cms.bool(True),
        select=cms.vstring('keep .* p4=12')
    ))
    process.acTriggerObjectSequence += cms.Sequence(getattr(process, 'ac_%s_ObjectProducer' % shortname))


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
        #'keep lheweights:lhapdf.306[0-9][0-9][0-9]=10',
        #'keep lheweights:PDF.306000=10',
        #'keep lheweights:dim6=10',
        #'keep lheweights:NNPDF31_nnlo_hessian_pdfas=10'
        )
)

process.acEventProducer = cms.EDProducer('AcornEventProducer')

process.p = cms.Path(
    process.egmPhotonIDSequence +
    process.selectedMuons +
    process.selectedPhotons +
    process.acMuonProducer +
    process.acPhotonProducer +
    process.acMCSequence +
    process.acTriggerObjectSequence +
    process.acEventInfoProducer +
    process.acEventProducer)

# process.schedule = cms.Schedule(process.patTriggerPath, process.p)
process.schedule = cms.Schedule(process.p)
# print process.dumpPython()
