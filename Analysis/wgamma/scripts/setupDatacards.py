#!/usr/bin/env python

import CombineHarvester.CombineTools.ch as ch
# import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--var', default='abs(reco_phi)')
parser.add_argument('--year', default='2016')
parser.add_argument('--label', default='default')
parser.add_argument('--channel', default='m')
parser.add_argument('--output', default='output/cards')
parser.add_argument('--type', default='eft_region', choices=['eft_region', 'pt_diff', 'fid_region'])
parser.add_argument('--pt-bins', type=int, default=6)
parser.add_argument('--phi-bins', type=int, default=5)
args = parser.parse_args()

cb = ch.CombineHarvester()

chn = args.channel
n_pt_bins = args.pt_bins
n_phi_bins = args.phi_bins
alt_shapes = True

charges = ['p', 'n']
fakes_name = 'data_fakes_highpt'
if args.type == 'fid_region':
    charges = ['x']
    fakes_name = 'data_fakes_sub'

era = '%s' % args.year
for c in charges:
    for ptbin in range(n_pt_bins):
        cat = (ptbin, '%s_%s_%i' % (c, chn, ptbin))
        cb.AddObservations(['*'], ['wg'], [era], [chn], [cat])
        cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_ooa_%s' % c, 'DY_XZG_R', 'ZG_IZG_R', 'DY_E', 'VV_R', 'VV_E', 'TT_XTTG_R', 'TTG_ITTG_R', 'TT_E', 'GG_R', 'GG_E', fakes_name], [cat], False)
        if args.type in ['eft_region', 'pt_phi_diff']:
            for phibin in range(n_phi_bins):
                cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_main_%s_%i_%i' % (c, ptbin, phibin)], [cat], True)
                cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_met1_%s_%i_%i' % (c, ptbin, phibin)], [cat], True)
        if args.type in ['pt_diff', 'fid_region']:
            for pt_truthbin in range(n_pt_bins):
                cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_main_%s_%i' % (c, pt_truthbin)], [cat], True)
                cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_met1_%s_%i' % (c, pt_truthbin)], [cat], True)

print '>> Adding systematic uncertainties...'


cb_no_fakes = cb.cp().process([fakes_name], False)

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
    for ix in xrange(1, 13):
        cb.cp().process(['data_fakes_sub']).AddSyst(
            cb, 'WeightStatSystBin%i' % ix, 'shape', ch.SystMap()(1.0))

##### TODO: QCD scale uncertainty
"""
    - In principle for all backgrounds (limited in practice)
    - For signal depends on measurement:
        * For 2D or 1D differential - only on met1 and OOA procs
        * For C3W - on everything
    - Should be de-correlated between procs
"""

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
cb.cp().ExtractShapes(
    file, '%s/$BIN/%s/$PROCESS' % (args.channel, args.var), '%s/$BIN/%s/$PROCESS_$SYSTEMATIC' % (args.channel, args.var))

cb.ForEachObj(lambda x: x.set_bin(x.bin() + '_' + args.year))


def DropBogusShapes(x):
    if x.type() == 'shape' and (x.value_u() <= 0. or x.value_d() <= 0.):
        print '>> Drop %s,%s,%s with val_u=%f, val_d=%f' % (x.bin(), x.process(), x.name(), x.value_u(), x.value_d())
        return True
    else:
        return False


cb.FilterSysts(lambda x: DropBogusShapes(x))


def CorrectionNegativeYield(proc):
    if proc.rate() < 0.:
        print '%s,%s has negative yield %f' % (proc.bin(), proc.process(), proc.rate())
        proc.set_rate(0.)


cb.ForEachProc(lambda x: CorrectionNegativeYield(x))
cb.SetAutoMCStats(cb, 0., True)
# cb.PrintAll()

writer = ch.CardWriter('$TAG/$BIN.txt',
                       '$TAG/common/$ANALYSIS.%s.%s.input.root' % (args.channel, args.year))
writer.SetWildcardMasses([])
writer.SetVerbosity(1)
writer.WriteCards(args.output, cb)

print '>> Done!'
