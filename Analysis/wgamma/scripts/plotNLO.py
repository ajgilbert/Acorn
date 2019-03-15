import ROOT
import argparse
import json
import itertools
from copy import deepcopy
from pprint import pprint
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
from array import array
from Acorn.Analysis.plottemplates import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
plot.ModTDRStyle(width=900, height=900)

hists = dict()
f_inc = ROOT.TFile('nlo_lhe_pt_inc.root')
hists['inc'] = f_inc.Get('lhe_pt')
f_120 = ROOT.TFile('nlo_lhe_pt_130_after_PS.root')
hists['120'] = f_120.Get('lhe_pt')

tot = hists['inc'].Integral(0, hists['inc'].GetNbinsX() + 1)
hi = hists['inc'].Integral(hists['inc'].GetXaxis().FindFixBin(120), hists['inc'].GetNbinsX() + 1)

print hi / tot

hists['inc'].Scale(1. / tot)
hists['120'].Scale((hi / tot) / hists['120'].Integral(0, hists['120'].GetNbinsX() + 1))

plotcfg = DEFAULT_CFG
plotcfg.update({
    'type': 'multihist',
    'legend_pos': [0.55, 0.66, 0.90, 0.91],
    'ratio_pad_frac': 0.33,
    'top_title_right': '41.5 fb^{-1} (13 TeV, 2017)',
    'logy': True,
    'logy_min': 1E-6
    })

MakeMultiHistPlot('nlo_lhe_pt_compare_130',
                    outdir='.',
                    hists=hists,
                    cfg=UpdateDict(plotcfg, {
                        'x_title': 'LHE photon p_{T} (GeV)',
                        'ratio_y_title': 'Ratio',
                        'top_title_right': '',
                        'ratio_y_range': [0.5, 1.5],
                        }),
                    layout=[
                        {'name': 'inc', 'legend': 'Inclusive'},
                        {'name': '120', 'legend': '120'},
                    ],
                    ratios=[
                        {'num': '120', 'den': 'inc'},
                    ])
# # Draw efficiencies for the stitched sample
# for year, var in itertools.product(years, ['l0_pt', 'p0_pt']):
#     hist_dict = {}
#     for opath, objname, obj in hists[year][var]['WG'].ListObjects(depth=0):
#         hist_dict[objname] = obj

#     # Draw N/N-1 efficiencies
#     MakeMultiHistPlot('efficiencies_m_%s_%s' % (year, var),
#                       outdir=outdir,
#                       hists=hist_dict,
#                       cfg=UpdateDict(plotcfg, {
#                           'x_title': x_titles[var],
#                           'ratio_y_title': '#varepsilon(N/N-1)',
#                           'top_title_right': '%.1f fb^{-1} (13 TeV, %s)' % (lumi_fb[year], year),
#                           'ratio_y_range': [0.55, 1.05],
#                           }),
#                       layout=[
#                           {'name': 'baseline_m', 'legend': 'Baseline'},
#                           {'name': 'muon_id', 'legend': ' + Muon ID'},
#                           {'name': 'muon_iso', 'legend': ' + Muon Iso'},
#                           {'name': 'muon_trg', 'legend': ' + Muon Trg(24/27) OR (50)'},
#                           {'name': 'muon_trg_2', 'legend': ' + Muon Trg(24/27)'},
#                           {'name': 'photon_id_m', 'legend': ' + Photon ID'},
#                           {'name': 'photon_pix_m', 'legend': ' + Photon pix veto'},
#                       ],
#                       ratios=[
#                           {'num': 'muon_id', 'den': 'baseline_m', 'type': 'binomial'},
#                           {'num': 'muon_iso', 'den': 'muon_id', 'type': 'binomial'},
#                           {'num': 'muon_trg', 'den': 'muon_iso', 'type': 'binomial'},
#                           {'num': 'muon_trg_2', 'den': 'muon_trg', 'type': 'binomial'},
#                           {'num': 'photon_id_m', 'den': 'muon_trg_2', 'type': 'binomial'},
#                           {'num': 'photon_pix_m', 'den': 'photon_id_m', 'type': 'binomial'},
#                       ]
#     )

