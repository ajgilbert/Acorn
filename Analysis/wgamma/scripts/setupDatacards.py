#!/usr/bin/env python

import CombineHarvester.CombineTools.ch as ch
# import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--var', default='abs(reco_phi)')
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

for c in ['p', 'n']:
    for ptbin in range(n_pt_bins):
        cat = (ptbin, '%s_m_%i' % (c, ptbin))
        cb.AddObservations(['*'], ['wg'], ["13TeV"], ['mu_%s' % c], [cat])
        cb.AddProcesses(['*'], ['wg'], ["13TeV"], ['mu_%s' % c], ['WG_%s_ooa' % c, 'VV_R', 'DY_R', 'TT_R', 'data_fakes'], [cat], False)
        for phibin in range(5):
            cb.AddProcesses(['*'], ['wg'], ["13TeV"], ['mu_%s' % c], ['WG_%s_%i_%i' % (c, ptbin, phibin)], [cat], True)
            cb.AddProcesses(['*'], ['wg'], ["13TeV"], ['mu_%s' % c], ['WG_%s_met1_%i_%i' % (c, ptbin, phibin)], [cat], True)

print '>> Adding systematic uncertainties...'

cb.cp().AddSyst(
    cb, 'lumi_$ERA', 'lnN', ch.SystMap()(1.026))

cb.cp().signals().AddSyst(
    cb, 'eff_$BIN', 'lnN', ch.SystMap()(1.1))

cb.cp().process(['WG_p_ooa', 'WG_n_ooa', 'DY_R', 'TTG_R', 'VV_R']).AddSyst(
    cb, 'eff_$BIN', 'lnN', ch.SystMap()(1.1))

cb.cp().process(['data_fakes']).AddSyst(
    cb, 'fakes_$BIN', 'lnN', ch.SystMap()(1.1))

cb.cp().AddSyst(
    cb, 'lumiscale', 'rateParam', ch.SystMap()(1.0))
par = cb.GetParameter('lumiscale').set_frozen(True)

print '>> Extracting histograms from input root files...'
file = 'output_2016_eft_region.root'
cb.cp().ExtractShapes(
    file, '$BIN/%s/$PROCESS' % args.var, '$BIN/%s/$PROCESS_$SYSTEMATIC' % args.var)

# cb.SetAutoMCStats(cb, 0, True)

writer = ch.CardWriter('$TAG/$ANALYSIS_$CHANNEL_$BINID_$ERA.txt',
                       '$TAG/common/$ANALYSIS.input.root')
writer.SetWildcardMasses([])
writer.SetVerbosity(1)
writer.WriteCards(args.output, cb)

print '>> Done!'