import FWCore.ParameterSet.Config as cms
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
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
opts.register('genOnly', 0, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Process as genOnly")
opts.register('updateJECs', 1, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Update JECs for jets and MET")
opts.register('hasLHE', 1, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Assume MC sample has LHE info")
opts.register('keepLHEWeights', 1, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Store some LHE weights")
opts.register('keepLHEParticles', 1, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Store the LHE particles")
opts.register('keepScaleWeights', 1, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Store the LHE scale weights")
opts.register('keepPDFWeights', 1, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Store the LHE pdf weights")
opts.register('doWGammaRivet', 0, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Run the WGamma RIVET routine and save output variables")
opts.register('cores', 1, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Number of cores/threads")
opts.register('input', 'root://xrootd.unl.edu//store/data/Run2016H/Tau/MINIAOD/PromptReco-v3/000/284/036/00000/36B9BD65-5B9F-E611-820B-02163E0126D3.root', parser.VarParsing.multiplicity.singleton, parser.VarParsing.varType.string, "input file")
opts.register('year', '2016', parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.string, "Year label")

opts.parseArguments()
isData = bool(opts.isData)
isMC = not isData
genOnly = int(opts.genOnly)
year = str(opts.year)
hasLHE = bool(opts.hasLHE)
keepLHEWeights = bool(opts.keepLHEWeights)
keepScaleWeights = bool(opts.keepScaleWeights)
keepPDFWeights = bool(opts.keepPDFWeights)
keepLHEParticles = bool(opts.keepLHEParticles)
updateJECs = bool(opts.updateJECs)
doWGammaRivet = bool(opts.doWGammaRivet)

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
# process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#         ignoreTotal = cms.untracked.int32(1)
#         )
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
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(opts.input.split(','))
)

process.GlobalTag.globaltag = cms.string(opts.globalTag)


process.options.numberOfThreads = cms.untracked.uint32(opts.cores)
process.options.numberOfStreams = cms.untracked.uint32(opts.cores)

process.selectedElectrons = cms.EDFilter("PATElectronRefSelector",
    src=cms.InputTag("slimmedElectrons"),
    cut=cms.string("pt > 10 & abs(eta) < 2.6")
)

process.selectedMuons = cms.EDFilter("PATMuonRefSelector",
    src=cms.InputTag("slimmedMuons"),
    cut=cms.string("pt > 9 & abs(eta) < 2.6")
)

process.selectedPhotons = cms.EDFilter("PATPhotonRefSelector",
    src=cms.InputTag("slimmedPhotons"),
    cut=cms.string("pt > 25 & abs(eta) < 2.6")
)

jetLabel = 'patJetsReapplyJECUpdatedJECs' if updateJECs else 'slimmedJets'
process.selectedJets = cms.EDFilter("PATJetRefSelector",
    src=cms.InputTag(jetLabel),
    cut=cms.string("pt > 30 & abs(eta) < 5.0")
)

process.acMuonProducer = cms.EDProducer('AcornMuonProducer',
    input=cms.InputTag("selectedMuons"),
    branch=cms.string('muons'),
    inputVertices=cms.InputTag('offlineSlimmedPrimaryVertices'),
    select=cms.vstring('keep .* p4=12 dxy=12 dz=12 vertex=12 pfIso.*=12')
)

process.customInitialSeq = cms.Sequence()

### Update the jets & MET in all years
if updateJECs:
    if year in ['2017']:
        # This is the EE noise mitigation - for 2017 data & MC only!
        # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_or_10
        runMetCorAndUncFromMiniAOD(
                process,
                isData=isData,
                fixEE2017=True,
                fixEE2017Params={'userawPt': True, 'ptThreshold': 50.0, 'minEtaThreshold': 2.65, 'maxEtaThreshold': 3.139},
                postfix="UpdatedJECs"
        )
    else:
        runMetCorAndUncFromMiniAOD(
                process,
                isData=isData,
                postfix="UpdatedJECs"
        )


    makePuppiesFromMiniAOD(process, True)
    runMetCorAndUncFromMiniAOD(process,
                               isData=isData,
                               metType="Puppi",
                               postfix="PuppiUpdatedJECs",
                               jetFlavor="AK4PFPuppi",
                               )
    process.puppiNoLep.useExistingWeights = True
    process.puppi.useExistingWeights = True
    # The sequence for PUPPI has to go first, then the normal MET, because of some overriding conflicts
    process.customInitialSeq += cms.Sequence(process.egmPhotonIDSequence + process.puppiMETSequence + process.fullPatMetSequencePuppiUpdatedJECs)
    process.customInitialSeq += cms.Sequence(process.fullPatMetSequenceUpdatedJECs)

### Setup for EGamma recipes
if year in ['2016_old']:
    setupEgammaPostRecoSeq(process, runEnergyCorrections=True, runVID=True, era='2016-Legacy')
elif year in ['2016']:
    setupEgammaPostRecoSeq(process, runEnergyCorrections=False, era='2016-Legacy')
elif year in ['2017']:
    setupEgammaPostRecoSeq(process, runVID=True, era='2017-Nov17ReReco')
elif year in ['2018']:
    setupEgammaPostRecoSeq(process, runEnergyCorrections=True, era='2018-Prompt')

### Electrons
ele_veto_id = "cutBasedElectronID-Fall17-94X-V2-veto"
ele_loose_id = "cutBasedElectronID-Fall17-94X-V2-loose"
ele_medium_id = "cutBasedElectronID-Fall17-94X-V2-medium"
ele_tight_id = "cutBasedElectronID-Fall17-94X-V2-tight"
ele_mva_wp80_id = "mvaEleID-Fall17-iso-V2-wp80"
ele_mva_wp90_id = "mvaEleID-Fall17-iso-V2-wp90"
ele_heep_id = "heepElectronID-HEEPV70"

# Embed the full vid results in the new electron object
process.slimmedElectrons.modifierConfig.modifications.append(cms.PSet(
    electron_config=cms.PSet(
        ElectronCutValues=cms.InputTag("egmGsfElectronIDs", ele_medium_id),
        electronSrc=(cms.InputTag("updatedElectrons") if year in ['2016'] else cms.InputTag("slimmedElectrons", "", "@skipCurrentProcess")),
    ),
    modifierName=cms.string('EGExtraInfoModifierFromVIDCutFlowResultValueMaps'),
    overrideExistingValues=cms.bool(True),
    photon_config=cms.PSet()
    )
)

for pset in process.egmGsfElectronIDs.physicsObjectIDs:
    idname = str(pset.idDefinition.idName.value())
    if 'cutBasedElectronID-Fall17-94X-V2-' in idname:
        for cut_pset in pset.idDefinition.cutFlow:
            if cut_pset.cutName.value() == 'GsfEleRelPFIsoScaledCut':
                print '>> Switching off isolation cut for %s' % idname
                cut_pset.isIgnored = cms.bool(True)

process.acElectronProducer = cms.EDProducer('AcornElectronProducer',
    input=cms.InputTag("selectedElectrons"),
    inputVertices=cms.InputTag('offlineSlimmedPrimaryVertices'),
    branch=cms.string('electrons'),
    select=cms.vstring('keep .* p4=12 dxy=12 dz=12 vertex=12 relativeEAIso=12 energyCorrections=12 scEta=12 scEnergy=12'),
    eleVetoIdMap=cms.InputTag(ele_veto_id),
    eleLooseIdMap=cms.InputTag(ele_loose_id),
    eleMediumIdMap=cms.InputTag(ele_medium_id),
    eleTightIdMap=cms.InputTag(ele_tight_id),
    eleMVAwp80IdMap=cms.InputTag(ele_mva_wp80_id),
    eleMVAwp90IdMap=cms.InputTag(ele_mva_wp90_id),
    eleHEEPIdMap=cms.InputTag(ele_heep_id),
    relativeEAIsoFromUserData=cms.vstring('ElectronCutValues', 'GsfEleRelPFIsoScaledCut_0'),
    takeIdsFromObjects=cms.bool(True),
    energyCorrections=cms.vstring('ecalTrkEnergyPostCorr', 'energyScaleUp', 'energyScaleDown', 'energySigmaUp', 'energySigmaDown')
)

### Photons
photon_loose_id = "cutBasedPhotonID-Fall17-94X-V2-loose"
photon_medium_id = "cutBasedPhotonID-Fall17-94X-V2-medium"
photon_tight_id = "cutBasedPhotonID-Fall17-94X-V2-tight"

process.acPhotonProducer = cms.EDProducer('AcornPhotonProducer',
    input=cms.InputTag("selectedPhotons"),
    branch=cms.string('photons'),
    select=cms.vstring('keep .* p4=12 scEta=12 hadTowOverEm=12 full5x5SigmaIetaIeta=12 .*Iso=12 energyCorrections=12'),
    phoLooseIdMap=cms.InputTag(photon_loose_id),
    phoMediumIdMap=cms.InputTag(photon_medium_id),
    phoTightIdMap=cms.InputTag(photon_tight_id),
    phoCutFlow=cms.InputTag(photon_medium_id),
    chargedIsolation=cms.string('phoChargedIsolation'),
    neutralHadronIsolation=cms.string('phoNeutralHadronIsolation'),
    photonIsolation=cms.string('phoPhotonIsolation'),
    worstChargedIsolation=cms.string('phoWorstChargedIsolation'),
    takeIdsFromObjects=cms.bool(True),
    energyCorrections=cms.vstring('ecalEnergyPostCorr', 'energyScaleUp', 'energyScaleDown', 'energySigmaUp', 'energySigmaDown')
)

### PFJets
process.acPFJetProducer = cms.EDProducer('AcornPFJetProducer',
    input=cms.InputTag("selectedJets"),
    branch=cms.string('pfJets'),
    select=cms.vstring('keep .* p4=12'),
    jetID=cms.PSet(
        version=cms.string('WINTER17'),
        quality=cms.string('TIGHT')
    )
)

### MET
metLabel = 'slimmedMETsUpdatedJECs' if updateJECs else 'slimmedMETs'
process.acPFType1MetProducer = cms.EDProducer('AcornMetProducer',
    input=cms.InputTag(metLabel),
    branch=cms.string('pfType1Met'),
    select=cms.vstring('keep .* p4=12', 'drop sumEt'),
    saveGenMetFromPat=cms.bool(False),
    saveCorrectionLevels=cms.vint32(1),
    saveUncertaintyShifts=cms.vint32(2, 3, 10, 11, 12, 13),
    skipMainMet=cms.bool(True)
)

puppiMetLabel = 'slimmedMETsPuppiUpdatedJECs' if updateJECs else 'slimmedMETsPuppi'
process.acPuppiMetProducer = cms.EDProducer('AcornMetProducer',
    input=cms.InputTag(puppiMetLabel),
    branch=cms.string('puppiMet'),
    select=cms.vstring('keep .* p4=12', 'drop sumEt'),
    saveGenMetFromPat=cms.bool(False),
    saveCorrectionLevels=cms.vint32(1),
    saveUncertaintyShifts=cms.vint32(2, 3, 10, 11, 12, 13),
    skipMainMet=cms.bool(True)
)

process.acMCSequence = cms.Sequence(
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.selectedGenParticles = cms.EDProducer("GenParticlePruner",
    src=cms.InputTag("prunedGenParticles"),
    select=cms.vstring(
        "drop  *",
        "keep isPromptFinalState()",
        "keep++ abs(pdgId) == 15 && statusFlags().isPrompt()",
        "+keep pdgId == 22 && status == 1 && (pt > 10 || isPromptFinalState())"
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
    select=cms.vstring('keep .* p4=12'),
    incomingP4Fix=cms.bool(True)
)

process.acPileupInfoProducer = cms.EDProducer('AcornPileupInfoProducer',
    input=cms.InputTag("slimmedAddPileupInfo"),
    branch=cms.string('pileupInfo'),
    select=cms.vstring('keep .*'),
    minBx=cms.int32(0),
    maxBx=cms.int32(0)
)

process.acGenMetProducer = cms.EDProducer('AcornMetProducer',
    input=cms.InputTag("genMetTrue"),
    branch=cms.string('genMet'),
    select=cms.vstring('keep .* p4=12'),
    saveGenMetFromPat=cms.bool(False),
    saveCorrectionLevels=cms.vint32(),
    saveUncertaintyShifts=cms.vint32(),
    skipMainMet=cms.bool(False)
)

process.acGenMetFromPatProducer = cms.EDProducer('AcornMetProducer',
    input=cms.InputTag("slimmedMETs"),
    branch=cms.string('genMet'),
    select=cms.vstring('keep .* p4=12'),
    saveGenMetFromPat=cms.bool(True),
    saveCorrectionLevels=cms.vint32(),
    saveUncertaintyShifts=cms.vint32(),
    skipMainMet=cms.bool(False)
)

process.selectedGenJets = cms.EDFilter("GenJetRefSelector",
    src=cms.InputTag("slimmedGenJets"),
    cut=cms.string("pt > 15")
)

process.acGenJetProducer = cms.EDProducer('AcornCandidateProducer',
    input=cms.InputTag("selectedGenJets"),
    branch=cms.string('genJets'),
    select=cms.vstring('keep .* p4=12')
)


if doWGammaRivet:
    process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
        inputPruned=cms.InputTag("prunedGenParticles"),
        inputPacked=cms.InputTag("packedGenParticles"),
    )
    process.generator = cms.EDProducer("GenParticles2HepMCConverter",
        genParticles=cms.InputTag("mergedGenParticles"),
        genEventInfo=cms.InputTag("generator", "", "SIM"),
        signalParticlePdgIds=cms.vint32()
    )
    process.wgammaRivetProducer = cms.EDProducer("WGammaRivetProducer",
        HepMCCollection=cms.InputTag("generator:unsmeared")
    )
    process.acWGammaRivetProducer = cms.EDProducer('AcornWGammaRivetProducer',
        input=cms.InputTag("wgammaRivetProducer"),
        branch=cms.string('rivetVariables'),
        select=cms.vstring('keep .*')
    )
    process.acMCSequence += cms.Sequence(
        process.mergedGenParticles +
        process.generator +
        process.wgammaRivetProducer +
        process.acWGammaRivetProducer
    )

if isMC:
    process.acMCSequence += cms.Sequence(
        process.selectedGenParticles +
        process.acGenParticleProducer +
        process.acGenMetFromPatProducer +
        process.acPileupInfoProducer +
        process.selectedGenJets +
        process.acGenJetProducer
    )
    if hasLHE and keepLHEParticles:
        process.acMCSequence += cms.Sequence(
            process.acLHEParticleProducer
        )

hlt_paths = [
    'HLT_IsoMu24_v',
    'HLT_IsoTkMu24_v',
    'HLT_IsoMu27_v',
    'HLT_Mu50_v',
    'HLT_TkMu50_v',
    'HLT_Ele27_WPTight_Gsf_v',
    'HLT_Ele32_WPTight_Gsf_L1DoubleEG_v',
    'HLT_Ele32_WPTight_Gsf_v',
    'HLT_Ele35_WPTight_Gsf_v',
    'HLT_Ele115_CaloIdVT_GsfTrkIdT_v',
    'HLT_Photon175_v',
    'HLT_Photon200_v',
]

process.acTriggerObjectSequence = cms.Sequence(
    # process.icTriggerPathProducer
)

trigger_objects = 'slimmedPatTrigger'
if year == '2016_old':
    trigger_objects = 'selectedPatTrigger'

for path in hlt_paths:
    shortname = path[4:-2]  # drop the HLT_ and _v parts
    setattr(process, 'ac_%s_ObjectProducer' % shortname, cms.EDProducer('AcornTriggerObjectProducer',
        input=cms.InputTag(trigger_objects),
        triggerResults=cms.InputTag('TriggerResults', '', 'HLT'),
        hltConfigProcess=cms.string('HLT'),
        branch=cms.string('triggerObjects_%s' % shortname),
        hltPath= cms.string(path),
        storeIfFired=cms.bool(True),
        select=cms.vstring('keep .* p4=12'),
        extraFilterLabels=cms.vstring()
    ))
    process.acTriggerObjectSequence += cms.Sequence(getattr(process, 'ac_%s_ObjectProducer' % shortname))

# Adds an extra filter label for the HLT_Ele32_WPTight_Gsf_L1DoubleEG path, needed to
# emulate HLT_Ele32_WPTight_Gsf in the portion of 2017 data that doesn't have it
# More detail: https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary#2017
# and slide 22 of: https://indico.cern.ch/event/662751/contributions/2778365/attachments/1561439/2458438/egamma_workshop_triggerTalk.pdf
process.ac_Ele32_WPTight_Gsf_L1DoubleEG_ObjectProducer.extraFilterLabels = cms.vstring('hltEGL1SingleEGOrFilter')


metfilter_proc_data = {
    "2016": "RECO",
    "2017": "PAT",
    "2018": "RECO"
}
metfilter_proc_mc = {
    "2016_old": "PAT",
    "2016": "PAT",
    "2017": "PAT",
    "2018": "PAT"
}
metfilter_proc = metfilter_proc_data if isData else metfilter_proc_mc


# Recommended extra metfilter to run on top of minioad in 2017 and 2018
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist,
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )


process.acEventInfoProducer = cms.EDProducer('AcornEventInfoProducer',
    lheProducer=cms.InputTag("externalLHEProducer"),
    generator=cms.InputTag("generator"),
    includeLHEInfo=cms.bool(isMC and hasLHE),
    includeLHEWeights=cms.bool(isMC and hasLHE and keepLHEWeights),
    includeGenWeights=cms.bool(isMC),
    metFilterResults=cms.InputTag("TriggerResults", "", metfilter_proc[year]),
    saveMetFilters=cms.vstring(
        'Flag_goodVertices',
        'Flag_globalSuperTightHalo2016Filter',
        'Flag_HBHENoiseFilter',
        'Flag_HBHENoiseIsoFilter',
        'Flag_EcalDeadCellTriggerPrimitiveFilter',
        'Flag_BadPFMuonFilter',
        'Flag_BadChargedCandidateFilter',
        'Flag_eeBadScFilter'
        ),
    saveMetFilterBools=cms.VInputTag(),
    userDoubles=cms.VInputTag(
        cms.InputTag('fixedGridRhoFastjetAll')
        ),
    branch=cms.string('eventInfo'),
    select=cms.vstring(
        'keep .* genWeights=10',
        'drop lheweights:.*',
        'drop lheweightgroups:.*',
        'keep lheweights:dim6'
        ),
    includeNumVertices=cms.bool(True),
    inputVertices=cms.InputTag('offlineSlimmedPrimaryVertices')
)


if keepScaleWeights:
    process.acEventInfoProducer.select += cms.vstring(
        'keep lheweights:(renscfact|facscfact|muR|muF|mur|muf|MUR|MUF).*=6'
    )


if keepPDFWeights:
    if year in ['2016', '2016_old']:
        process.acEventInfoProducer.select += cms.vstring(
            'keep lheweightgroups:.*NNPDF30_lo_as_0130.LHgrid.*=6',  # 1
            'keep lheweightgroups:.*NNPDF30_lo_as_0130_nf_4.LHgrid.*=6',  # 2
            'keep lheweightgroups:.*NNPDF30_nlo_nf_5_pdfas.*=6',  # 3
            'keep lheweightgroups:.*NNPDF31_nnlo_as_0118.*=6',  # 4
            'keep lheweights:.*PDF.set...260[0-9][0-9][0-9].*=6',  # 5
        )
    if year in ['2017', '2018']:
        process.acEventInfoProducer.select += cms.vstring(
            # 'keep lheweightgroups:.*NNPDF31_nlo_as_0118_nf_4.*=6',  # 1
            'keep lheweightgroups:.*NNPDF31_nnlo_as_0118_nf_4.*=6',  # 2
            # 'keep lheweightgroups:.*NNPDF31_nlo_hessian_pdfas.*=6',  # 3
            'keep lheweightgroups:.*NNPDF31_nnlo_hessian_pdfas.*=6',  # 4
            'keep lheweights:.*lhapdf.306[0-9][0-9][0-9].*=6',  # 5
            # 'keep lheweights:.*lhapdf.305[0-9][0-9][0-9].*=6',  # 6
            'keep lheweightgroups:.*NNPDF31_nnlo_as_0118_mc_hessian_pdfas.*=6',  # 7
        )


if year in ['2017', '2018'] and genOnly == 0:
    process.customInitialSeq += process.ecalBadCalibReducedMINIAODFilter
    process.acEventInfoProducer.saveMetFilterBools = cms.VInputTag(
        cms.InputTag('ecalBadCalibReducedMINIAODFilter')
        )


### Prefiring weights - only needed for 2016 and 2017 MC
prefiring_era_str = {
    "2016_old": "2016BtoH",
    "2016": "2016BtoH",
    "2017": "2017BtoF",
    "2018": "2017BtoF"
}

from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
    DataEra = cms.string(prefiring_era_str[year]), #Use 2016BtoH for 2016
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUncty = cms.double(0.2),
    SkipWarnings = False
)

if isMC and year in ['2016_old', '2016', '2017'] and genOnly == 0:
    process.customInitialSeq += cms.Sequence(process.prefiringweight)
    process.acEventInfoProducer.userDoubles += cms.VInputTag(
        cms.InputTag('prefiringweight:nonPrefiringProb'),
        cms.InputTag('prefiringweight:nonPrefiringProbUp'),
        cms.InputTag('prefiringweight:nonPrefiringProbDown')
        )

process.acEventProducer = cms.EDProducer('AcornEventProducer')

if genOnly == 1:
    # Take the full collection for now
    process.acGenParticleProducer.input = cms.InputTag("prunedGenParticles")
    process.acEventInfoProducer.saveMetFilters = cms.vstring()
    process.acEventInfoProducer.userDoubles = cms.VInputTag()
    process.selectedGenJets.src = cms.InputTag("ak4GenJetsNoNu")
    process.acEventInfoProducer.includeNumVertices=cms.bool(False)
    process.p = cms.Path(
        process.acGenMetProducer +
        process.acGenParticleProducer +
        process.acLHEParticleProducer +
        process.selectedGenJets +
        process.acGenJetProducer +
        process.acEventInfoProducer +
        process.acEventProducer)
elif genOnly == 2:
    # Take the full collection for now
    process.acEventInfoProducer.saveMetFilters = cms.vstring()
    process.acEventInfoProducer.includeNumVertices=cms.bool(False)
    process.acEventInfoProducer.userDoubles = cms.VInputTag()
    process.p = cms.Path(
        process.acLHEParticleProducer +
        process.acEventInfoProducer +
        process.acEventProducer)
else:
    process.p = cms.Path(
        process.egammaPostRecoSeq +
        process.customInitialSeq +
        process.selectedElectrons +
        process.selectedMuons +
        process.selectedPhotons +
        process.selectedJets +
        process.acElectronProducer +
        process.acMuonProducer +
        process.acPhotonProducer +
        process.acPFJetProducer +
        process.acPFType1MetProducer +
        process.acPuppiMetProducer +
        process.acMCSequence +
        process.acTriggerObjectSequence +
        process.acEventInfoProducer +
        process.acEventProducer)

# process.schedule = cms.Schedule(process.patTriggerPath, process.p)
process.schedule = cms.Schedule(process.p)
# print process.dumpPython()
