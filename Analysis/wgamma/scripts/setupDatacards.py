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
wg_proc = 'WG_NLO'
alt_shapes = False

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

# https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#SummaryTable
cb.cp().AddSyst(
    cb, 'lumi_$ERA', 'lnN', ch.SystMap('era')
    (['2016'], 1.023)
    (['2017'], 1.025)
    (['2018'], 1.023))

cb_no_fakes = cb.cp().process([fakes_name], False)

cb_no_fakes.cp().AddSyst(
    cb, 'CMS_eff_$CHANNEL', 'shape', ch.SystMap()(1.0))

cb_no_fakes.cp().AddSyst(
    cb, 'CMS_trigger_$CHANNEL', 'shape', ch.SystMap()(1.0))

cb_no_fakes.cp().AddSyst(
    cb, 'CMS_eff_p', 'shape', ch.SystMap()(1.0))

cb_no_fakes.cp().AddSyst(
    cb, 'CMS_ele_fake_p', 'shape', ch.SystMap()(1.0))

cb_no_fakes.cp().AddSyst(
    cb, 'CMS_prefiring', 'shape', ch.SystMap()(1.0))

if fakes_name == 'data_fakes_sub':
    for ix in xrange(1, 13):
        cb.cp().process(['data_fakes_sub']).AddSyst(
            cb, 'data_fakes_WeightStatSystBin%i' % ix, 'shape', ch.SystMap()(1.0))

if alt_shapes:
    cb_no_fakes.cp().AddSyst(
        cb, 'CMS_scale_p', 'shape', ch.SystMap()(1.0))

    cb_no_fakes.cp().AddSyst(
        cb, 'CMS_scale_met_jes', 'shape', ch.SystMap()(1.0))

    cb_no_fakes.cp().AddSyst(
        cb, 'CMS_scale_met_unclustered', 'shape', ch.SystMap()(1.0))


cb.cp().AddSyst(
    cb, 'lumiscale', 'rateParam', ch.SystMap()(1.0))
par = cb.GetParameter('lumiscale').set_frozen(True)

print '>> Extracting histograms from input root files...'
file = 'output_%s_%s_%s.root' % (args.year, args.type, args.label)
cb.cp().ExtractShapes(
    file, '%s/$BIN/%s/$PROCESS' % (args.channel, args.var), '%s/$BIN/%s/$PROCESS_$SYSTEMATIC' % (args.channel, args.var))

# cb.SetAutoMCStats(cb, 0, True)
cb.ForEachObj(lambda x: x.set_bin(x.bin() + '_' + args.year))

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
