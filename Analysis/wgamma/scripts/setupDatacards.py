#!/usr/bin/env python

import CombineHarvester.CombineTools.ch as ch
# import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--var', default='abs(reco_phi)')
parser.add_argument('--year', default='2016')
parser.add_argument('--channel', default='m')
parser.add_argument('--output', default='output/cards')
args = parser.parse_args()

cb = ch.CombineHarvester()


# procs = {
#   'sig'  : ['tH_YtMinus', 'tHW'],
#   'sim'  : ['WZ', 'ZZ', 'ttW', 'ttZ', 'ttH'],
#   'bkg'  : ['WZ','ZZ','ttW','ttZ','ttH','reducible']
# }

# cats = [(0, 'emt'), (1, 'mmt')]

# masses = ['125']
n_pt_bins = 6
era = '13TeV_%s' % args.year
for c in ['p', 'n']:
    for ptbin in range(n_pt_bins):
        cat = (ptbin, '%s_%s_%i' % (c, args.channel, ptbin))
        cb.AddObservations(['*'], ['wg'], [era], ['mu_%s' % c], [cat])
        cb.AddProcesses(['*'], ['wg'], [era], ['mu_%s' % c], ['WG_%s_ooa' % c, 'DY_R', 'VV_R', 'TTG_R', 'data_fakes'], [cat], False)
        for phibin in range(5):
            cb.AddProcesses(['*'], ['wg'], [era], ['mu_%s' % c], ['WG_%s_%i_%i' % (c, ptbin, phibin)], [cat], True)
            cb.AddProcesses(['*'], ['wg'], [era], ['mu_%s' % c], ['WG_%s_met1_%i_%i' % (c, ptbin, phibin)], [cat], True)

print '>> Adding systematic uncertainties...'

cb.cp().AddSyst(
    cb, 'lumi_$ERA', 'lnN', ch.SystMap()(1.026))

cb.cp().signals().AddSyst(
    cb, 'eff_$BIN_$ERA', 'lnN', ch.SystMap()(1.1))

cb.cp().process(['WG_p_ooa', 'WG_n_ooa', 'VV_R', 'DY_R', 'TTG_R']).AddSyst(
    cb, 'eff_$BIN_$ERA', 'lnN', ch.SystMap()(1.1))

cb.cp().process(['data_fakes']).AddSyst(
    cb, 'fakes_$BIN_$ERA', 'lnN', ch.SystMap()(1.1))

cb.cp().AddSyst(
    cb, 'lumiscale', 'rateParam', ch.SystMap()(1.0))
par = cb.GetParameter('lumiscale').set_frozen(True)

print '>> Extracting histograms from input root files...'
file = 'output_%s_eft_region.root' % args.year
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
