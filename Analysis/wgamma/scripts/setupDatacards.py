#!/usr/bin/env python

import CombineHarvester.CombineTools.ch as ch
# import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--var', default='abs(reco_phi)')
parser.add_argument('--year', default='2016')
parser.add_argument('--channel', default='m')
parser.add_argument('--output', default='output/cards')
parser.add_argument('--type', default='eft', choices=['eft', 'pt_diff'])
args = parser.parse_args()

cb = ch.CombineHarvester()


# procs = {
#   'sig'  : ['tH_YtMinus', 'tHW'],
#   'sim'  : ['WZ', 'ZZ', 'ttW', 'ttZ', 'ttH'],
#   'bkg'  : ['WZ','ZZ','ttW','ttZ','ttH','reducible']
# }

# cats = [(0, 'emt'), (1, 'mmt')]

# masses = ['125']
chn = args.channel
n_pt_bins = 6
wg_proc = 'WG_NLO'

era = '13TeV_%s' % args.year
for c in ['p', 'n']:
    for ptbin in range(n_pt_bins):
        cat = (ptbin, '%s_%s_%i' % (c, chn, ptbin))
        cb.AddObservations(['*'], ['wg'], [era], [chn], [cat])
        cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_%s_ooa' % c, 'DY_XZG_R', 'ZG_IZG_R', 'DY_E', 'VV_R', 'VV_E', 'TT_R', 'TT_E', 'GG', 'data_fakes'], [cat], False)
        if args.type in ['eft']:
            for phibin in range(5):
                cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_%s_%i_%i' % (c, ptbin, phibin)], [cat], True)
                cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_%s_met1_%i_%i' % (c, ptbin, phibin)], [cat], True)
        if args.type in ['pt_diff']:
            for pt_truthbin in range(n_pt_bins):
                cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_%s_met1_%i' % (c, pt_truthbin)], [cat], True)
                cb.AddProcesses(['*'], ['wg'], [era], [chn], ['WG_%s_%i' % (c, pt_truthbin)], [cat], True)

print '>> Adding systematic uncertainties...'

cb.cp().AddSyst(
    cb, 'lumi_$ERA', 'lnN', ch.SystMap()(1.026))

cb.cp().signals().AddSyst(
    cb, 'eff_$ERA', 'lnN', ch.SystMap()(1.1))

cb.cp().process(['WG_p_ooa', 'WG_n_ooa', 'VV_R', 'DY_XZG_R', 'ZG_IZG_R', 'TT_R']).AddSyst(
    cb, 'eff_$ERA', 'lnN', ch.SystMap()(1.1))

cb.cp().process(['data_fakes']).AddSyst(
    cb, 'fakes_$ERA', 'lnN', ch.SystMap()(1.1))

cb.cp().AddSyst(
    cb, 'lumiscale', 'rateParam', ch.SystMap()(1.0))
par = cb.GetParameter('lumiscale').set_frozen(True)

print '>> Extracting histograms from input root files...'
file = 'output_%s_eft_region_default.root' % args.year
cb.cp().ExtractShapes(
    file, '%s/$BIN/%s/$PROCESS' % (args.channel, args.var), '%s/$BIN/%s/$PROCESS_$SYSTEMATIC' % (args.channel, args.var))

# cb.SetAutoMCStats(cb, 0, True)
cb.ForEachObj(lambda x: x.set_bin(x.bin() + '_' + args.year))

writer = ch.CardWriter('$TAG/$BIN.txt',
                       '$TAG/common/$ANALYSIS.%s.%s.input.root' % (args.channel, args.year))
writer.SetWildcardMasses([])
writer.SetVerbosity(1)
writer.WriteCards(args.output, cb)

print '>> Done!'
