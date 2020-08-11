#!/usr/bin/env python
import ROOT
import CombineHarvester.CombineTools.ch as ch
# import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--var', default='abs(reco_phi)')
parser.add_argument('--year', default='2016')
parser.add_argument('--label', default='default')
parser.add_argument('--channel', default='m')
parser.add_argument('--output', default='output/cards')
parser.add_argument('--type', default='eft_region', choices=['baseline', 'eft_region', 'pt_diff', 'fid_region'])
parser.add_argument('--pt-bins', type=int, default=6)
parser.add_argument('--phi-bins', type=int, default=5)
args = parser.parse_args()

cb = ch.CombineHarvester()
cb.SetFlag('filters-use-regex', True)

chn = args.channel
n_pt_bins = args.pt_bins
n_phi_bins = args.phi_bins
alt_shapes = True

N_FAKE_BINS = 18
DO_WORST_ISO = 2  # 0 = No, 1 = Only the cut for photon fakes, 2 = full
if DO_WORST_ISO >= 2:
    N_FAKE_BINS = 16

charges = ['x']
if args.type == 'eft_region' and 'chg' in args.label:
    charges = ['p', 'n']
fakes_name = 'data_fakes_highpt'
if args.type in ['baseline', 'fid_region']:
    charges = ['x']
    fakes_name = 'data_fakes_sub'

era = '%s' % args.year
for c in charges:
    for ptbin in range(n_pt_bins):
        cat = (ptbin, '%s_%s_%i' % (c, chn, ptbin))
        cb.AddObservations(['*'], ['wg'], [era], [chn], [cat])
        cb.AddProcesses(['*'], ['wg'], [era], [chn], ['DY_XZG_R', 'ZG_IZG_R', 'DY_E', 'VV_R', 'VV_E', 'ST_R', 'ST_E', 'TT_XTTG_R', 'TTG_ITTG_R', 'TT_E', 'GG_R', 'GG_E', fakes_name], [cat], False)
        if args.type in ['baseline']:
            cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_R'], [cat], True)
        if args.type in ['fid_region', 'eft_region']:
            cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_ooa_%s' % c], [cat], False)
        if args.type in ['baseline', 'fid_region']:
            cb.AddProcesses(['*'], ['wg'], [era], [chn], ['data_fakes_lep_sub'], [cat], False)
        if args.type in ['eft_region', 'pt_phi_diff']:
            if ptbin == 0:
                pt_truthbins = [ptbin, ptbin + 1]
            elif ptbin == n_pt_bins - 1:
                pt_truthbins = [ptbin - 1, ptbin]
            else:
                pt_truthbins = [ptbin - 1, ptbin, ptbin + 1]
            for pt_truthbin in pt_truthbins:
                for phibin in range(n_phi_bins):
                    cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_main_%s_%i_%i' % (c, pt_truthbin, phibin)], [cat], True)
                    cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_met1_%s_%i_%i' % (c, pt_truthbin, phibin)], [cat], True)
        if args.type in ['pt_diff', 'fid_region']:
            for pt_truthbin in range(n_pt_bins):
                cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_main_%s_%i' % (c, pt_truthbin)], [cat], True)
                cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_met1_%s_%i' % (c, pt_truthbin)], [cat], True)

print '>> Adding systematic uncertainties...'


cb_no_fakes = cb.cp().process([fakes_name, 'data_fakes_lep_sub'], False)

##### Lumi
# https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#SummaryTable
# Uncorrelated 2016,2.2,0.0,0.0
# Uncorrelated 2017,0.0,2.0,0.0
# Uncorrelated 2018,0.0,0.0,1.5
# X-Y factorization,0.9,0.8,2.0
# Length scale 17-18,0.0,0.3,0.2
# Beam-beam deflection 15-17,0.4,0.4,0.0
# Dynamic beta 15-17,0.5,0.5,0.0
# Beam current calibration 17-18,0.0,0.3,0.2
# Ghosts and satellites 15-17,0.4,0.1,0.0
cb_no_fakes.cp().AddSyst(
    cb, 'CMS_lumi_$ERA', 'lnN', ch.SystMap('era')
    (['2016'], 1.022)
    (['2017'], 1.020)
    (['2018'], 1.015))

cb_no_fakes.cp().AddSyst(
    cb, 'CMS_lumi_xy_fact', 'lnN', ch.SystMap('era')
    (['2016'], 1.009)
    (['2017'], 1.008)
    (['2018'], 1.020))

cb_no_fakes.cp().AddSyst(
    cb, 'CMS_lumi_length_scale', 'lnN', ch.SystMap('era')
    (['2017'], 1.003)
    (['2018'], 1.002))

cb_no_fakes.cp().AddSyst(
    cb, 'CMS_lumi_deflection', 'lnN', ch.SystMap('era')
    (['2016'], 1.004)
    (['2017'], 1.004))

cb_no_fakes.cp().AddSyst(
    cb, 'CMS_lumi_dynamic_beta', 'lnN', ch.SystMap('era')
    (['2016'], 1.005)
    (['2017'], 1.005))

cb_no_fakes.cp().AddSyst(
    cb, 'CMS_lumi_beam_current', 'lnN', ch.SystMap('era')
    (['2017'], 1.003)
    (['2018'], 1.002))

cb_no_fakes.cp().AddSyst(
    cb, 'CMS_lumi_ghosts', 'lnN', ch.SystMap('era')
    (['2016'], 1.004)
    (['2017'], 1.001))

##### TODO: Pileup uncertainty
cb_no_fakes.cp().AddSyst(
    cb, 'CMS_pileup', 'shape', ch.SystMap()(1.0))

##### Electron & muon ID & iso efficiencies
cb_no_fakes.cp().AddSyst(
    cb, 'CMS_eff_$CHANNEL', 'shape', ch.SystMap()(1.0))

##### Electron & muon trigger efficiences
cb_no_fakes.cp().AddSyst(
    cb, 'CMS_trigger_$CHANNEL', 'shape', ch.SystMap()(1.0))

##### Photon ID & iso efficiency
cb_no_fakes.cp().AddSyst(
    cb, 'CMS_eff_p', 'shape', ch.SystMap()(1.0))

##### Electron -> photon fake rate
cb_no_fakes.cp().AddSyst(
    cb, 'CMS_ele_fake_p', 'shape', ch.SystMap()(1.0))

##### Prefiring
cb_no_fakes.cp().AddSyst(
    cb, 'CMS_prefiring', 'shape', ch.SystMap()(1.0))

##### Photon fake estimation (low pT method only)
if fakes_name == 'data_fakes_sub':
    for ix in xrange(1, N_FAKE_BINS + 1):
        cb.cp().process(['data_fakes_sub']).AddSyst(
            cb, 'WeightStatSystBin%i' % ix, 'shape', ch.SystMap()(1.0))

##### Lepton fake estimation
cb.cp().process(['data_fakes_lep_sub']).AddSyst(
    cb, 'CMS_fake_$CHANNEL_$ERA_norm', 'lnN', ch.SystMap('channel', 'era')
    (['e'], ['2016'], 1.44)
    (['e'], ['2017'], 1.25)
    (['e'], ['2018'], 1.31)
    (['m'], ['2016'], 1.29)
    (['m'], ['2017'], 1.23)
    (['m'], ['2018'], 1.18)
    )

if args.type in ['baseline']:
    cb.cp().process(['WG_R']).AddSyst(
        cb, 'QCDScale', 'shape', ch.SystMap()(1.0))

##### QCD scale uncertainty
if args.type in ['fid_region']:
    cb.cp().process(['WG_main_.*']).AddSyst(
        cb, 'QCDScaleAccept', 'shape', ch.SystMap()(1.0))
    cb.cp().process(['WG_met1_.*']).AddSyst(
        cb, 'QCDScaleRatio', 'shape', ch.SystMap()(1.0))
    cb.cp().process(['WG_ooa_.*']).AddSyst(
        cb, 'QCDScale', 'shape', ch.SystMap()(1.0))

if args.type in ['eft_region']:
    # Only adding scale uncerts on the diagonal for now...
    for ptbin in range(n_pt_bins):
        cb.cp().bin(['._._%i' % ptbin]).process(['WG_main_._%i_.*' % ptbin]).AddSyst(
            cb, 'QCDScale', 'shape', ch.SystMap()(1.0))
        cb.cp().bin(['._._%i' % ptbin]).process(['WG_met1_._%i_.*' % ptbin]).AddSyst(
            cb, 'QCDScaleRatio', 'shape', ch.SystMap()(1.0))
    cb.cp().process(['WG_ooa_.*']).AddSyst(
        cb, 'QCDScale', 'shape', ch.SystMap()(1.0))

"""
    - In principle for all backgrounds (limited in practice)
    - For signal depends on measurement:
        * For 2D or 1D differential - only on met1 and OOA procs
        * For C3W - on everything
    - Should be de-correlated between procs
"""

cb.cp().process(['DY_XZG_R', 'DY_E']).AddSyst(
    cb, 'QCDScale_Z_NNLO', 'lnN', ch.SystMap()(1.02))

cb.cp().process(['ZG_IZG_R']).AddSyst(
    cb, 'QCDScale_ZG_NLO', 'lnN', ch.SystMap()(1.035))

cb.cp().process(['VV_R', 'VV_E']).AddSyst(
    cb, 'QCDScale_VV_NLO', 'lnN', ch.SystMap()(1.034))

cb.cp().process(['TT_XTTG_R', 'TT_E']).AddSyst(
    cb, 'QCDScale_ttbar_NNLO', 'lnN', ch.SystMap()(1.034))

cb.cp().process(['GG_R', 'GG_E']).AddSyst(
    cb, 'QCDScale_GG_NLO', 'lnN', ch.SystMap()(1.16))

cb.cp().process(['VV_R', 'VV_E']).AddSyst(
    cb, 'pdf_VV', 'lnN', ch.SystMap()(1.048))

cb.cp().process(['TT_XTTG_R', 'TT_E']).AddSyst(
    cb, 'pdf_ttbar', 'lnN', ch.SystMap()(1.042))

# WW: 3% scale, 4.6% pdf
# WZ: 3.4% scale, 4.7% pdf
# ZZ: 3% scale, 4.8% pdf
# ST t-channel: 3% scale, 3% pdf
# ST tW-channel: 2.5% scale, 4.3% pdf
# ST s-channel: 3% scale, 3% pdf
# WWG?
# TGJets?

# ZG: 3.5% from scale variations

if alt_shapes:
    ##### Photon energy scale
    cb_no_fakes.cp().AddSyst(
        cb, 'CMS_scale_p', 'shape', ch.SystMap()(1.0))

    ##### MET scale due to JES
    cb_no_fakes.cp().AddSyst(
        cb, 'CMS_scale_met_jes', 'shape', ch.SystMap()(1.0))

    ##### MET scale due to unclustered energy
    cb_no_fakes.cp().AddSyst(
        cb, 'CMS_scale_met_unclustered', 'shape', ch.SystMap()(1.0))

    ##### TODO: electron energy scale

    ##### TODO: muon momentum scale

cb.cp().AddSyst(
    cb, 'lumiscale', 'rateParam', ch.SystMap()(1.0))
par = cb.GetParameter('lumiscale').set_frozen(True)

print '>> Extracting histograms from input root files...'
file = 'output_%s_%s_%s_merged.root' % (args.year, args.type, args.label)
if not alt_shapes:
    file = 'output_%s_%s_%s.root' % (args.year, args.type, args.label)

if args.type == 'baseline':
    cb.cp().ExtractShapes(
        file, '%s/baseline_%s_mZ_veto/%s/$PROCESS' % (args.channel, args.channel, args.var), '%s/baseline_%s_mZ_veto/%s/$PROCESS_$SYSTEMATIC' % (args.channel, args.channel, args.var))
else:
    cb.cp().ExtractShapes(
        file, '%s/$BIN/%s/$PROCESS' % (args.channel, args.var), '%s/$BIN/%s/$PROCESS_$SYSTEMATIC' % (args.channel, args.var))

cb.ForEachObj(lambda x: x.set_bin(x.bin() + '_' + args.year))

decorrelate_years = [
    'CMS_eff_e',
    'CMS_eff_m',
    'CMS_trigger_e',
    'CMS_trigger_m',
    'CMS_eff_p',
    'CMS_ele_fake_p',
    'WeightStatSystBin.*',
    'CMS_scale_p',
    'CMS_scale_met_jes',
    'CMS_scale_met_unclustered'
]

if args.type in ['fid_region', 'baseline']:
    decorrelate_proc = [
        'QCDScale',
        'QCDScaleAccept',
        'QCDScaleRatio'
    ]
    decorrelate_pt_bin = []
elif args.type in ['eft_region']:
    decorrelate_proc = [
        'QCDScale'
    ]

    decorrelate_pt_bin = [
        'QCDScale',
        'QCDScaleAccept',
        'QCDScaleRatio',
    ]

for syst_pattern in decorrelate_years:
    for syst_name in cb.cp().syst_name([syst_pattern]).syst_name_set():
        for era in cb.cp().syst_name([syst_name]).SetFromSysts(lambda x: x.era()):
            cb.cp().era([era]).RenameSystematic(cb, syst_name, syst_name + '_' + era)
            cb.RenameParameter(syst_name, syst_name + '_' + era)

for syst_pattern in decorrelate_pt_bin:
    for syst_name in cb.cp().signals().syst_name([syst_pattern]).syst_name_set():
        for proc in cb.cp().signals().syst_name([syst_name]).SetFromSysts(lambda x: x.process()):
            pt_bin = proc.split('_')[3]
            newname = syst_name + '_WG_PtBin_%s' % pt_bin
            cb.cp().process([proc]).RenameSystematic(cb, syst_name, newname)
            cb.RenameParameter(syst_name, newname)

### Fully decorrelates based on the process name, but perhaps overkill?
for syst_pattern in decorrelate_proc:
    for syst_name in cb.cp().syst_name([syst_pattern]).syst_name_set():
        for proc in cb.cp().syst_name([syst_name]).SetFromSysts(lambda x: x.process()):
            cb.cp().process([proc]).RenameSystematic(cb, syst_name, syst_name + '_' + proc)
            cb.RenameParameter(syst_name, syst_name + '_' + proc)


# cb.cp().syst_name(decorrelate_years).ForEachSyst(lambda x: x.set_name(x.name() + '_' + x.era()))
# cb.cp().syst_name(decorrelate_proc).ForEachSyst(lambda x: x.set_name(x.name() + '_' + x.process()))


def DropBogusShapes(x):
    if x.type() == 'shape' and (x.value_u() <= 0. or x.value_d() <= 0.):
        print '>> Drop %s,%s,%s with val_u=%f, val_d=%f' % (x.bin(), x.process(), x.name(), x.value_u(), x.value_d())
        return True
    else:
        return False


cb.FilterSysts(lambda x: DropBogusShapes(x))

# Can safely do this if every channel is a 1-bin counting expt.
if args.type in ['fid_region']:
    cb.cp().syst_type(['shape']).ForEachSyst(lambda x: x.set_type('lnN'))


def CorrectionNegativeYield(proc):
    if proc.rate() < 0.:
        print '%s,%s has negative yield %f' % (proc.bin(), proc.process(), proc.rate())
        proc.set_rate(0.)


def CollectZeroYields(proc, res):
    if proc.rate() == 0.:
        res.add((proc.bin(), proc.process()))


def SetDummyTemplate(proc):
    if proc.rate() == 0.:
        h = proc.ShapeAsTH1F()
        h.SetBinContent(int(h.GetNbinsX() / 2), 1E-7)
        h.Print('range')
        proc.set_shape(h, True)


cb.ForEachProc(lambda x: CorrectionNegativeYield(x))

zero_yields = set()
cb.ForEachProc(lambda x: CollectZeroYields(x, zero_yields))
print zero_yields
if args.type in 'baseline':
    cb.ForEachProc(lambda x: SetDummyTemplate(x))
    cb.FilterSysts(lambda x: (x.bin(), x.process()) in zero_yields)
    cb.ForEachObj(lambda x: x.set_bin('%s_%s_%s' % (args.channel, args.var, args.year)))
else:
    cb.FilterSysts(lambda x: (x.bin(), x.process()) in zero_yields)
    cb.FilterProcs(lambda x: (x.bin(), x.process()) in zero_yields)


cb.SetAutoMCStats(cb, 0., True)
# cb.PrintAll()

# Groups:
# - stat, lumi, theory, expt
# - stat, lumi, theory, mcstats, lepton_eff, photon_eff, photon_scale, photon_fakes, met_scale
cb.SetGroup('all', ['.*'])

gr_lumi = ['CMS_lumi_.*']
gr_lep_eff = ['CMS_eff_[em]_.*', 'CMS_trigger_[em]_.*']
gr_pho_eff = ['CMS_eff_p_.*']
gr_pho_fakes = ['CMS_ele_fake_p_.*', 'WeightStatSystBin.*']
gr_lep_fakes = ['CMS_fake_[em]_.*']
gr_other = ['CMS_pileup', 'CMS_prefiring']
gr_met_scale = ['CMS_scale_met_.*']
gr_pho_scale = ['CMS_scale_p_.*']
gr_sig_th = ['QCDScale.*_WG_.*']
gr_sig_inc = ['QCDScale_WG_PtBin_.*']
gir_bkg_th = ['QCDScale_.*NLO', 'pdf_VV', 'pdf_ttbar']

gr_th = gr_sig_th + gir_bkg_th
gr_expt = gr_lep_eff + gr_pho_eff + gr_pho_fakes + gr_lep_fakes + gr_other + gr_met_scale + gr_pho_scale

cb.SetGroup('lumi', gr_lumi)
cb.SetGroup('lep_eff', gr_lep_eff)
cb.SetGroup('pho_eff', gr_pho_eff)
cb.SetGroup('pho_fakes', gr_pho_fakes)
cb.SetGroup('lep_fakes', gr_lep_fakes)
cb.SetGroup('other', gr_other)
cb.SetGroup('met_scale', gr_met_scale)
cb.SetGroup('pho_scale', gr_pho_scale)
cb.SetGroup('sig_th', gr_sig_th)
cb.SetGroup('sig_inc', gr_sig_inc)
cb.SetGroup('bkg_th', gir_bkg_th)
cb.SetGroup('th', gr_th)
cb.SetGroup('expt', gr_expt)

cb.PrintParams()

writer = ch.CardWriter('$TAG/$BIN.txt',
                       '$TAG/common/$ANALYSIS.%s.%s.input.root' % (args.channel, args.year))
writer.SetWildcardMasses([])
writer.SetVerbosity(1)
writer.WriteCards(args.output, cb)

print '>> Done!'
